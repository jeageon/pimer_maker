# Primer Maker

`pimer_maker`는 기존 `DH5a-UTG` 기반 코드 구조를 활용해
GenBank(`.gb`) 파일을 입력받아
`프라이머 간섭 구간`과 `프라이머 후보`를 계산하고
결과를 주석(Feature) 포함 `GenBank`로 내려주는 도구입니다.

## 실행 방법

### macOS / Linux

```bash
cd /Users/jg/Documents/prime_maker/pimer_maker
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m streamlit run src/webui.py
```

또는:

```bash
./scripts/launch_primer_webui.command
```

### Windows (더블클릭 실행)

1. `scripts\launch_primer_webui.bat` 실행
2. Python 3.10+가 설치되어 있는지 확인
3. 의존성 설치 후 브라우저에서 Streamlit 화면이 열립니다.

#### Windows 바탕화면 바로가기 만들기

```bash
scripts\make_primer_desktop_shortcut.bat
```

`make_primer_desktop_shortcut.bat`를 실행하면 바탕화면에 `Primer Maker.lnk`가 생성됩니다.

## 사용 흐름

1. GenBank 업로드
2. 템플릿 간섭 영역 추출 (기존 feature + extreme GC / homopolymer / ambiguous)
3. 18~24 nt 프라이머를 Tm/GC/3' clamp/반복서열 기준으로 스캔
4. 통과 후보를 feature로 `GenBank`에 기록
5. 결과 파일 다운로드
6. 최종 페어는 사용자가 직접 고른 뒤 하단 입력으로 간섭 확인

## 핵심 파일

- `src/webui.py`: 업로드/파이프라인 실행 UI
- `src/modules/primer_pipeline.py`: 파이프라인 핵심 로직
  - 간섭 영역 수집
  - 후보 프라이머 탐색
  - 후보 프라이머 간섭 점검
  - GB 주석 생성
- `scripts/launch_primer_webui.bat`: Windows 더블클릭 실행기
- `scripts/make_primer_desktop_shortcut.bat`: 바탕화면 바로가기 생성기
- `requirements.txt`: 기존 `biopython`, `streamlit` 기반 의존성

## GitHub 새 저장소 업데이트 가이드

1. GitHub에서 `pimer_maker` 저장소 생성
2. 로컬에서 프로젝트 파일을 push
3. `.bat` 실행 스크립트와 `README.md` 포함 여부 확인

예:

```bash
git init
git add .
git commit -m "Init primer maker"
git remote add origin <YOUR_GITHUB_REPO_URL>
git branch -M main
git push -u origin main
```
