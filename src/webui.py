from __future__ import annotations

from src.config import (
    DEFAULT_GC_CLAMP_MAX,
    DEFAULT_GC_CLAMP_MIN,
    DEFAULT_GC_MAX,
    DEFAULT_GC_MIN,
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
st.caption("최종 쌍 선택은 사람 검토 후 아래에서 확인만 수행")

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

if st.button("쌍 검증"):
    if product_min >= product_max:
        st.error("Amplicon 최소 크기는 최대 크기보다 작아야 합니다.")
        st.stop()
    if not forward_input.strip() or not reverse_input.strip():
        st.error("두 프라이머 시퀀스를 모두 입력하세요.")
        st.stop()
    check = validate_primer_pair(
        sequence=result["sequence"],
        forward_seq=forward_input,
        reverse_seq=reverse_input,
        product_min=int(product_min),
        product_max=int(product_max),
        tm_target=float(tm_target),
        tm_tolerance=float(tm_tolerance),
    )
    if not check["valid"]:
        st.error("검증 실패")
        for msg in check.get("errors", []):
            st.warning(msg)
    else:
        st.success("조건 통과")
        if check.get("warnings"):
            for msg in check["warnings"]:
                st.warning(msg)
        st.dataframe(check["pairs"], use_container_width=True, hide_index=True)

with st.expander("실행 메타데이터"):
    st.json(result["metadata"])
