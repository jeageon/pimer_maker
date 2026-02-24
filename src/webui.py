from __future__ import annotations

from typing import Any

from src.config import (
    DEFAULT_GC_CLAMP_MAX,
    DEFAULT_GC_CLAMP_MIN,
    DEFAULT_GC_MAX,
    DEFAULT_GC_MIN,
    DEFAULT_MANUAL_HAIRPIN_MAX_K,
    DEFAULT_MANUAL_HAIRPIN_MIN_K,
    DEFAULT_MANUAL_OFFTARGET_SEED_LEN,
    DEFAULT_MANUAL_OFFTARGET_SEED_WARNING_LIMIT,
    DEFAULT_MANUAL_PAIR_DIMER_MAX_OVERLAP,
    DEFAULT_MANUAL_PAIR_DIMER_MIN_OVERLAP,
    DEFAULT_MANUAL_REQUIRE_3P_DIMER,
    DEFAULT_MANUAL_SELF_DIMER_MAX_OVERLAP,
    DEFAULT_MANUAL_SELF_DIMER_MIN_OVERLAP,
    DEFAULT_MANUAL_TM_GAP_FAIL,
    DEFAULT_MAX_PRODUCT_SIZE,
    DEFAULT_MIN_PRODUCT_SIZE,
    DEFAULT_PRIMER_LEN_MAX,
    DEFAULT_PRIMER_LEN_MIN,
    DEFAULT_PRIMER_TM_TARGET,
    DEFAULT_PRIMER_TM_TOLERANCE,
)
from src.modules.primer_pipeline import run_primer_pipeline, validate_primer_pair
import streamlit as st


st.set_page_config(page_title="Primer Maker", page_icon="🧬", layout="wide")
st.title("Primer Maker")
st.caption("GenBank (.gb) 업로드 기반 프라이머 전체 후보 설계")


uploaded_file = st.file_uploader("GenBank 파일", type=["gb", "gbk", "genbank"])
previous_result = st.session_state.get("last_result")
result = None


def _normalize_primer_result(raw_result: object) -> dict[str, object] | None:
    if not isinstance(raw_result, dict):
        return None
    metadata = raw_result.get("metadata")
    if not isinstance(metadata, dict):
        return None
    required_keys = ("metadata", "sequence", "gb_text", "filename")
    for key in required_keys:
        if key not in raw_result:
            return None
    return raw_result


def _dual_lang_text(message: str) -> str:
    translations: dict[str, str] = {
        "both forward and reverse primer sequences are required": "정방향/역방향 프라이머 시퀀스가 모두 필요합니다.",
        "primers must contain only A/C/G/T/N characters": "프라이머에 허용되지 않는 문자가 포함되어 있습니다. A/C/G/T/N만 허용됩니다.",
        "failed to calculate Tm": "Tm 계산에 실패했습니다.",
        "Tm gap too large": "정방향/역방향 Tm 차이가 허용치보다 큽니다.",
        "forward primer is likely to form hairpin": "정방향 프라이머가 hairpin(자체 접힘) 위험이 높습니다.",
        "reverse primer is likely to form hairpin": "역방향 프라이머가 hairpin(자체 접힘) 위험이 높습니다.",
        "forward primer is likely to form self-dimer": "정방향 프라이머가 self-dimer(자체 이량체) 위험이 높습니다.",
        "reverse primer is likely to form self-dimer": "역방향 프라이머가 self-dimer(자체 이량체) 위험이 높습니다.",
        "forward primer does not match template": "정방향 프라이머가 템플릿에서 매칭되는 위치가 없습니다.",
        "reverse primer does not match template": "역방향 프라이머가 템플릿에서 매칭되는 위치가 없습니다.",
        "no valid product found under configured range": "현재 설정 조건에서 조건을 만족하는 증폭 산물이 없습니다.",
        "all opposite-strand pairs were outside product range": "반대 가닥 조합의 모든 쌍이 Amplicon 크기 조건에서 제외되었습니다.",
        "some opposite-strand pairs were filtered by product size range": "반대 가닥 조합의 일부가 Amplicon 크기 조건에 의해 제외되었습니다.",
        "all opposite-strand pairs were rejected by dimer risk filter": "반대 가닥 조합의 모든 쌍이 이량체 위험 필터에서 제외되었습니다.",
        "some opposite-strand pairs were rejected by dimer risk filter": "반대 가닥 조합의 일부가 이량체 위험 필터에서 제외되었습니다.",
        "all candidate pairs show high off-target seed-seed potential": "모든 후보 조합에서 오프타겟 seed 위험도가 높게 계산되었습니다.",
        "candidate pair(s) show high off-target seed-seed potential": "일부 후보 조합에서 오프타겟 seed 위험도가 높습니다.",
        "no valid pair found after current validation filters": "현재 검증 필터 조건을 모두 통과하는 쌍이 없습니다.",
        "forward/reverse binding sites are on the same strand; reverse primer may need reverse-complement input": "정방향/역방향 결합 위치가 동일 가닥입니다. 역방향 프라이머는 상보서열(reverse complement)로 입력이 필요할 수 있습니다.",
        "forward_tm_outside_target_": "정방향 Tm이 목표 범위를 벗어났습니다.",
        "reverse_tm_outside_target_": "역방향 Tm이 목표 범위를 벗어났습니다.",
        "Tm gap": "Tm 차이",
        "exceeds tolerance": "허용 Tm 오차를 초과했습니다.",
        "outside target": "목표 Tm 범위를 벗어났습니다.",
        "no candidate pairs passed quality filters": "모든 후보 쌍이 간섭/오차 기준을 통과하지 못했습니다.",
        "pair": "쌍(Forward/Reverse)",
        "pair binding positions are too short": "쿼리 간 접촉이 짧아 신뢰성이 낮습니다.",
        "forward primer has no perfect match in template": "정방향 프라이머가 템플릿에서 완전 일치 위치를 찾지 못했습니다.",
        "reverse primer has no perfect match in template": "역방향 프라이머가 템플릿에서 완전 일치 위치를 찾지 못했습니다.",
    }
    for key, value in translations.items():
        if message == key or message.startswith(key):
            return f"{message}\n- 한국어: {value}"
    return message


def _extract_filter_summary_rows(filter_summary: object) -> list[dict[str, str]]:
    if not isinstance(filter_summary, dict):
        return []
    return [
        {"구분": "총 후보 결합쌍", "값": str(filter_summary.get("total_combinations", 0)), "비고": "forward x reverse 모든 결합조합"},
        {"구분": "반대가닥 결합쌍", "값": str(filter_summary.get("opposite_pairs_total", 0)), "비고": "정방향/역방향 가닥이 반대인 조합"},
        {"구분": "same-strand 제외", "값": str(filter_summary.get("skipped_same_strand", 0)), "비고": "같은 가닥만 매칭된 조합"},
        {"구분": "Amplicon 크기 제외", "값": str(filter_summary.get("skipped_size", 0)), "비고": "설정한 Amplicon 범위 밖"},
        {"구분": "이량체 필터 제외", "값": str(filter_summary.get("skipped_dimer", 0)), "비고": "cross-dimer/hairpin 기준 제외"},
        {"구분": "오프타겟 경고 표시", "값": str(filter_summary.get("off_target_pairs", 0)), "비고": f"seed_len={filter_summary.get('seed_len', '-')}, limit={filter_summary.get('warning_limit', '-')}"}
    ]


def _render_filter_status(filter_summary: dict[str, object]) -> None:
    if not filter_summary:
        return
    total = int(filter_summary.get("total_combinations", 0) or 0)
    opposite = int(filter_summary.get("opposite_pairs_total", 0) or 0)
    skipped_size = int(filter_summary.get("skipped_size", 0) or 0)
    skipped_dimer = int(filter_summary.get("skipped_dimer", 0) or 0)
    off_target = int(filter_summary.get("off_target_pairs", 0) or 0)

    if opposite <= 0:
        st.error("반대가닥 결합쌍이 없습니다.")
    if skipped_size == opposite and opposite > 0:
        st.error("Amplicon 크기 조건에서 모든 반대가닥 결합쌍이 제외되었습니다.")
    elif skipped_size > 0:
        st.warning(f"Amplicon 크기 조건에서 {skipped_size}개가 제외되었습니다.")

    if skipped_dimer == opposite and opposite > 0:
        st.error("이량체 위험 필터에서 반대가닥 결합쌍이 모두 제외되었습니다.")
    elif skipped_dimer > 0:
        st.warning(f"이량체 위험 필터에서 {skipped_dimer}개가 제외되었습니다.")

    if off_target and opposite > 0:
        if off_target >= opposite:
            st.error("오프타겟 경고가 전체 결합쌍에서 높은 수준으로 감지되었습니다.")
        else:
            st.warning(f"오프타겟 경고가 {off_target}개 결합쌍에서 감지되었습니다.")

    if total > 0 and opposite <= 0 and skipped_size == 0 and skipped_dimer == 0:
        st.warning("same-strand 결합만 존재해 유효한 증폭 쌍을 만들지 못했습니다.")


def _render_interference_details(interference_details: object) -> None:
    if not isinstance(interference_details, list) or not interference_details:
        return

    rows: list[dict[str, str]] = []

    for item in interference_details:
        if not isinstance(item, dict):
            continue

        scope = str(item.get("scope", ""))
        inter_type = str(item.get("type", ""))
        primer_type = str(item.get("primer", ""))

        if scope == "single_primer":
            target = "정방향" if primer_type == "forward" else "역방향"
            if inter_type == "hairpin":
                rows.append(
                    {
                        "대상": target,
                        "간섭 유형": "자기접힘(헤어핀)",
                        "설명": f"{target} 말단이 내부 서열과 상보적으로 접힘 (길이 {item.get('overlap_len')}nt)",
                        "프라이머 말단": str(item.get("primer_tail", "")),
                        "내부 매칭 서열": str(item.get("internal_match", "")),
                        "내부 위치(0-base)": str(item.get("match_index", "")),
                    }
                )
            elif inter_type == "self_dimer":
                rows.append(
                    {
                        "대상": target,
                        "간섭 유형": "자기이량체",
                        "설명": f"같은 프라이머끼리 보체 결합 가능성 존재 (겹침 {item.get('overlap_len')}nt)",
                        "프라이머 말단": str(item.get("primer_tail", "")),
                        "상보 매칭": str(item.get("complement_match", "")),
                        "매칭 위치(0-base)": str(item.get("paired_index", "")),
                    }
                )
            else:
                rows.append({"대상": target, "간섭 유형": str(inter_type), "설명": json_like(item)})

        elif scope == "pair":
            if inter_type == "cross_dimer":
                rows.append(
                    {
                        "대상": "쌍(Forward/Reverse)",
                        "간섭 유형": "이량체(cross-dimer)",
                        "설명": f"Forward/Reverse 사이 상보 결합 가능 (겹침 {item.get('overlap_len')}nt)",
                        "Forward 위치(0-base)": str(item.get("forward_pos", "")),
                        "Reverse 위치(0-base)": str(item.get("reverse_pos", "")),
                        "Forward 말단": str(item.get("a_tail", "")),
                        "Reverse 말단": str(item.get("b_tail", "")),
                    }
                )
            elif inter_type == "off_target":
                rows.append(
                    {
                        "대상": "쌍(Forward/Reverse)",
                        "간섭 유형": "오프타겟 시드",
                        "설명": f"seed_len={item.get('seed_len', '')}, Forward/Reverse seed_hits={item.get('forward_seed_hits', '')}/{item.get('reverse_seed_hits', '')}",
                        "Forward seed": str(item.get("forward_seed", "")),
                        "Reverse seed": str(item.get("reverse_seed", "")),
                    }
                )
            elif inter_type == "tm_gap":
                rows.append(
                    {
                        "대상": "쌍(Forward/Reverse)",
                        "간섭 유형": "Tm 차이 과다",
                        "설명": f"ΔTm={float(item.get('gap', 0.0)):.2f}°C (허용치 {item.get('limit', 0.0)}°C)",
                        "forward Tm": str(item.get("forward_tm", "")),
                        "reverse Tm": str(item.get("reverse_tm", "")),
                    }
                )
            elif inter_type == "binding_miss":
                rows.append(
                    {
                        "대상": "쌍(Forward/Reverse)",
                        "간섭 유형": "결합 실패",
                        "설명": f"{'정방향' if primer_type == 'forward' else '역방향'} 프라이머가 템플릿에서 매칭되지 않음",
                    }
                )
            else:
                rows.append({"대상": "쌍(Forward/Reverse)", "간섭 유형": str(inter_type), "설명": json_like(item)})
        else:
            rows.append({"대상": "-", "간섭 유형": str(inter_type or scope), "설명": json_like(item)})

    if not rows:
        return

    st.subheader("간섭 시퀀스 상세")
    st.dataframe(rows[:20], use_container_width=True, hide_index=True)


def json_like(value: Any) -> str:
    if isinstance(value, dict):
        return ", ".join(f"{k}={v}" for k, v in value.items())
    return str(value)

with st.form("design_form"):
    st.subheader("설계 조건")
    c1, c2 = st.columns(2)
    with c1:
        tm_target = st.number_input("목표 Tm (°C)", 40.0, 80.0, DEFAULT_PRIMER_TM_TARGET, 0.5)
    with c2:
        tm_tolerance = st.number_input("Tm 허용 오차 (±°C)", 0.0, 15.0, DEFAULT_PRIMER_TM_TOLERANCE, 0.5)

    c4, c5, c6 = st.columns(3)
    with c4:
        len_min = st.number_input("최소 프라이머 길이", 14, 32, DEFAULT_PRIMER_LEN_MIN, 1)
    with c5:
        len_max = st.number_input("최대 프라이머 길이", 14, 32, DEFAULT_PRIMER_LEN_MAX, 1)
    with c6:
        repeat_run_limit = st.number_input("동일염기 반복 제한(권장 4)", 2, 10, 4, 1)

    c7, c8, c9 = st.columns(3)
    with c7:
        gc_min = st.number_input("GC min (%)", 0.0, 100.0, DEFAULT_GC_MIN, 1.0)
    with c8:
        gc_max = st.number_input("GC max (%)", 0.0, 100.0, DEFAULT_GC_MAX, 1.0)

    c10, c11 = st.columns(2)
    with c10:
        clamp_min = st.number_input("3' GC clamp min", 0, 5, DEFAULT_GC_CLAMP_MIN, 1)
    with c11:
        clamp_max = st.number_input("3' GC clamp max", 0, 5, DEFAULT_GC_CLAMP_MAX, 1)

    include_input_features = st.checkbox(
        "입력 GB feature를 간섭영역(회피 대상)으로 사용",
        value=False,
        help="체크하면 CDS/promoter 등 기존 어노테이션 구간에서 프라이머 후보를 제외합니다.",
    )

    run_btn = st.form_submit_button("프라이머 설계 실행")


if not run_btn and previous_result is None:
    if not uploaded_file:
        st.info("GenBank 파일을 업로드하고 실행 버튼을 눌러 주세요.")
    st.stop()

if run_btn and not uploaded_file:
    st.error("GenBank 파일이 필요합니다.")
    st.stop()

if run_btn and gc_min > gc_max:
    st.error("GC min은 GC max보다 작아야 합니다.")
    st.stop()
if run_btn and len_min > len_max:
    st.error("최소 길이는 최대 길이보다 작아야 합니다.")
    st.stop()

if run_btn:
    if gc_min > gc_max or len_min > len_max:
        st.stop()
    with st.spinner("간섭영역 분석 및 후보 생성 중..."):
        try:
            result = run_primer_pipeline(
                uploaded_file.getvalue(),
                source_filename=uploaded_file.name,
                tm_target=float(tm_target),
                tm_tolerance=float(tm_tolerance),
                gc_min=float(gc_min),
                gc_max=float(gc_max),
                len_min=int(len_min),
                len_max=int(len_max),
                repeat_run_limit=int(repeat_run_limit),
                gc_clamp_min=int(clamp_min),
                gc_clamp_max=int(clamp_max),
                include_input_features=bool(include_input_features),
            )
        except Exception as exc:
            st.error(f"실행 실패: {exc}")
            st.exception(exc)
            st.stop()
        if not isinstance(result, dict) or not isinstance(result.get("metadata"), dict):
            st.error("설계 결과가 올바른 형식이 아닙니다.")
            st.stop()
        st.session_state["last_result"] = result
else:
    normalized = _normalize_primer_result(previous_result)
    if normalized is not None:
        result = normalized
    else:
        result = None
if result is None:
    st.error("설계 결과가 없습니다. 먼저 '프라이머 설계 실행'을 눌러 결과를 생성해 주세요.")
    st.stop()

st.success("설계 완료")
result = _normalize_primer_result(result)
if result is None:
    st.error("저장된 결과 포맷이 손상되었습니다. 다시 설계를 실행해 주세요.")
    st.stop()
meta = result["metadata"]
lcol, rcol = st.columns(2)
with lcol:
    st.subheader("요약")
    display_name = result.get("record_name") or result.get("record_id")
    st.write(f"Record Name: `{display_name}`")
    st.write(f"시퀀스 길이: `{meta['sequence_length']}`")
    st.write(f"간섭영역 수: `{meta['interference_count']}`")
    st.write(f"프라이머 후보 수: `{meta['primer_candidate_count']}`")
with rcol:
    st.subheader("진행 파라미터")
    st.write(f"Tm: `{meta['tm_target']} ± {meta['tm_tolerance']}`")
    st.write(f"GC: `{meta['gc_min']}~{meta['gc_max']}`")
    st.write(f"길이: `{meta['len_min']}~{meta['len_max']}`")
    st.write(f"반복 제한: `{meta['repeat_run_limit']}`")
    st.write(f"간섭 구분: `{meta.get('interference_by_type', {})}`")

st.download_button(
    label="주석 포함 GenBank 다운로드",
    data=result["gb_text"],
    file_name=result["filename"],
    mime="application/octet-stream",
    use_container_width=True,
)

st.subheader("프라이머 후보 미리보기")
primer_candidates = result.get("primer_candidates")
if primer_candidates:
    st.dataframe(primer_candidates[:50], use_container_width=True, hide_index=True)
else:
    st.warning("조건을 만족하는 후보가 없습니다. 조건을 완화해 보세요.")


st.divider()
st.subheader("프라이머 쌍 간섭 확인 (수동)")
st.caption("최종 쌍 선택은 사용자 검토 후 아래에서 수용")

fc1, fc2 = st.columns(2)
with fc1:
    forward_input = st.text_area("Forward primer (5'-3')", height=80)
with fc2:
    reverse_input = st.text_area("Reverse primer (5'-3')", height=80)

pc1, pc2 = st.columns(2)
with pc1:
    product_min = st.number_input("Amplicon 최소", 30, 20000, DEFAULT_MIN_PRODUCT_SIZE, 10)
with pc2:
    product_max = st.number_input("Amplicon 최대", 50, 50000, DEFAULT_MAX_PRODUCT_SIZE, 10)

with st.expander("쌍 간섭 검사 고급 설정 (수동 검증)"):
    hc1, hc2 = st.columns(2)
    with hc1:
        pair_tm_gap_fail = st.number_input(
            "Tm 차이 실패 임계값 (°C)",
            0.0,
            20.0,
            DEFAULT_MANUAL_TM_GAP_FAIL,
            0.5,
            help="정방향/역방향 Tm 차이가 이 값보다 크면 즉시 실패 처리",
        )
    with hc2:
        pair_3p_anchor_only = st.checkbox(
            "교차 이량체 3' 말단 기준 필터링",
            value=DEFAULT_MANUAL_REQUIRE_3P_DIMER,
        )

    hc3, hc4 = st.columns(2)
    with hc3:
        pair_hairpin_min_k = st.number_input("Hairpin min k", 1, 12, DEFAULT_MANUAL_HAIRPIN_MIN_K, 1)
        pair_hairpin_max_k = st.number_input("Hairpin max k", 1, 12, DEFAULT_MANUAL_HAIRPIN_MAX_K, 1)
    with hc4:
        pair_self_dimer_min_overlap = st.number_input(
            "Self-dimer overlap min",
            1,
            12,
            DEFAULT_MANUAL_SELF_DIMER_MIN_OVERLAP,
            1,
        )
        pair_self_dimer_max_overlap = st.number_input(
            "Self-dimer overlap max",
            1,
            12,
            DEFAULT_MANUAL_SELF_DIMER_MAX_OVERLAP,
            1,
        )

    hc5, hc6 = st.columns(2)
    with hc5:
        pair_cross_dimer_min_overlap = st.number_input(
            "Cross-dimer overlap min",
            1,
            12,
            DEFAULT_MANUAL_PAIR_DIMER_MIN_OVERLAP,
            1,
        )
    with hc6:
        pair_cross_dimer_max_overlap = st.number_input(
            "Cross-dimer overlap max",
            1,
            12,
            DEFAULT_MANUAL_PAIR_DIMER_MAX_OVERLAP,
            1,
        )

    offt1, offt2 = st.columns(2)
    with offt1:
        pair_offtarget_seed_len = st.number_input(
            "오프타겟 seed 길이",
            4,
            12,
            DEFAULT_MANUAL_OFFTARGET_SEED_LEN,
            1,
        )
    with offt2:
        pair_offtarget_seed_warning_limit = st.number_input(
            "오프타겟 warning 임계값",
            0,
            20,
            DEFAULT_MANUAL_OFFTARGET_SEED_WARNING_LIMIT,
            1,
            help="0이면 warning 카운팅을 비활성화합니다.",
        )

if st.button("쌍 검증"):
    if product_min >= product_max:
        st.error("Amplicon 최소 크기는 최대 크기보다 작아야 합니다.")
        st.stop()
    if not forward_input.strip() or not reverse_input.strip():
        st.error("두 프라이머 시퀀스를 모두 입력하세요.")
        st.stop()
    if pair_hairpin_min_k > pair_hairpin_max_k:
        st.error("Hairpin min k는 max k보다 클 수 없습니다.")
        st.stop()
    if pair_self_dimer_min_overlap > pair_self_dimer_max_overlap:
        st.error("Self-dimer overlap min은 max보다 클 수 없습니다.")
        st.stop()
    if pair_cross_dimer_min_overlap > pair_cross_dimer_max_overlap:
        st.error("Cross-dimer overlap min은 max보다 클 수 없습니다.")
        st.stop()
    if pair_offtarget_seed_len < 4:
        st.error("오프타겟 seed 길이는 4 이상으로 입력해야 합니다.")
        st.stop()
    check = validate_primer_pair(
        sequence=result["sequence"],
        forward_seq=forward_input,
        reverse_seq=reverse_input,
        product_min=int(product_min),
        product_max=int(product_max),
        tm_gap_fail=float(pair_tm_gap_fail),
        hairpin_min_k=int(pair_hairpin_min_k),
        hairpin_max_k=int(pair_hairpin_max_k),
        self_dimer_min_overlap=int(pair_self_dimer_min_overlap),
        self_dimer_max_overlap=int(pair_self_dimer_max_overlap),
        pair_dimer_min_overlap=int(pair_cross_dimer_min_overlap),
        pair_dimer_max_overlap=int(pair_cross_dimer_max_overlap),
        pair_dimer_require_3p=bool(pair_3p_anchor_only),
        offtarget_seed_len=int(pair_offtarget_seed_len),
        offtarget_seed_warning_limit=int(pair_offtarget_seed_warning_limit),
        tm_target=float(tm_target),
        tm_tolerance=float(tm_tolerance),
    )
    if not check["valid"]:
        st.error("검증 실패")
        _render_interference_details(check.get("interference_details"))
        with st.expander("검증 결과 상세", expanded=True):
            for msg in check.get("summary_messages", []):
                translated = _dual_lang_text(str(msg))
                if "오류" in translated:
                    st.error(translated)
                elif "제외" in str(msg) or "표시" in str(msg) or "임계값" in str(msg):
                    st.warning(translated)
                else:
                    st.info(translated)
            filter_summary = check.get("filter_summary")
            if filter_summary:
                st.subheader("필터 요약")
                _render_filter_status(filter_summary)
                filter_rows = _extract_filter_summary_rows(filter_summary)
                if filter_rows:
                    st.dataframe(filter_rows, use_container_width=True, hide_index=True)
                size_examples = filter_summary.get("product_size_samples")
                if isinstance(size_examples, list) and size_examples:
                    st.caption(f"산물 크기 제외 예시: {size_examples[:10]}")
                off_target_examples = filter_summary.get("off_target_samples")
                if isinstance(off_target_examples, list) and off_target_examples:
                    st.caption(f"오프타겟 시드 예시(앞->뒤): {off_target_examples[:10]}")
                st.json(filter_summary)
        for msg in check.get("errors", []):
            st.error(_dual_lang_text(str(msg)))
    else:
        st.success("조건 통과")
        _render_interference_details(check.get("interference_details"))
        with st.expander("검증 결과 상세", expanded=False):
            for msg in check.get("summary_messages", []):
                translated = _dual_lang_text(str(msg))
                if "제외" in str(msg) or "경고" in str(msg) or "오류" in str(msg):
                    st.warning(translated)
                else:
                    st.info(translated)
            filter_summary = check.get("filter_summary")
            if filter_summary:
                st.subheader("필터 요약")
                _render_filter_status(filter_summary)
                filter_rows = _extract_filter_summary_rows(filter_summary)
                if filter_rows:
                    st.dataframe(filter_rows, use_container_width=True, hide_index=True)
                size_examples = filter_summary.get("product_size_samples")
                if isinstance(size_examples, list) and size_examples:
                    st.caption(f"산물 크기 제외 예시: {size_examples[:10]}")
                off_target_examples = filter_summary.get("off_target_samples")
                if isinstance(off_target_examples, list) and off_target_examples:
                    st.caption(f"오프타겟 시드 예시(앞->뒤): {off_target_examples[:10]}")
                st.json(filter_summary)
        if check.get("warnings"):
            for msg in check.get("warnings"):
                translated = _dual_lang_text(str(msg))
                if "Tm" in str(msg):
                    st.warning(translated)
                else:
                    st.info(translated)
        st.dataframe(check["pairs"], use_container_width=True, hide_index=True)

with st.expander("실행 메타데이터"):
    st.json(result["metadata"])
