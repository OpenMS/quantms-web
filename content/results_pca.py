"""PCA Results Page."""
import streamlit as st
import pandas as pd
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from src.common.common import page_setup
from src.common.results_helpers import get_abundance_data
from openms_insight.components.pca import run_and_plot_pca

params = page_setup()
st.title("PCA Analysis")

st.markdown(
    """
Principal Component Analysis (PCA) of protein-level abundance.
Samples are colored by group assignment to visualize clustering.
"""
)

if "workspace" not in st.session_state:
    st.warning("Please initialize your workspace first.")
    st.stop()

# 1. 변경된 get_abundance_data 적용 (반환값 2개: pivot_df, group_map)
result = get_abundance_data(st.session_state["workspace"])
if result is None:
    st.info("Abundance data not available. Please run the workflow and configure sample groups first.")
    st.page_link("content/results_abundance.py", label="Go to Abundance", icon="📋")
    st.stop()

_, group_map = result

if "statistics_df" not in st.session_state or st.session_state["statistics_df"] is None:
    st.info("Statistical analysis data not found. Please run the statistical inference first to obtain p-adj values.")
    # st.page_link("content/results_statistical.py", label="Go to Statistical Inference", icon="📊")
    st.stop()

target_df = st.session_state["statistics_df"]

# 2. 이 페이지에서 직접 expr_df(발현량 매트릭스) 구축하기
# group_map의 키(샘플명들)를 컬럼으로 사용하여 발현량 데이터만 추출합니다.
sample_columns = list(group_map.keys())

# pivot_df에 단백질 식별자(예: ProteinName)와 샘플 컬럼들이 포함되어 있어야 합니다.
if "ProteinName" in target_df.columns:
    expr_df = target_df.set_index("ProteinName")[sample_columns]
elif target_df.index.name == "ProteinName":
    expr_df = target_df[sample_columns]
else:
    # 예외 방지: ProteinName이 컬럼에 없고 인덱스 이름도 지정되지 않은 경우 첫 번째 컬럼을 인덱스로 가정
    expr_df = target_df.set_index(target_df.columns[0])[sample_columns]

top_n = 500

# 3. p-value 기준 상위 n개 단백질 필터링
# pivot_df에서 유의미한 단백질 탐색
top_proteins = (
    target_df
    .dropna(subset=["p-adj"])
    .sort_values("p-adj", ascending=True)
    .head(top_n)
)

# 만약 위에서 인덱스를 바꿨다면 pivot_df 구조에 맞게 단백질 이름을 가져옵니다.
if "ProteinName" in top_proteins.columns:
    top_protein_names = top_proteins["ProteinName"]
else:
    top_protein_names = top_proteins.index

expr_df_pca = expr_df.loc[
    expr_df.index.intersection(top_protein_names)
]

if expr_df_pca.shape[0] < 2:
    st.info("Not enough proteins after p-value filtering for PCA.")
    st.stop()

# 4. OpenMS-Insight 모듈 호출 파트
# 이 아래의 지저분한 계산 및 Plotly 시각화 코드를 외부 모듈로 캡슐화하여 호출합니다.
try:
    # 정의된 분석 및 시각화 함수 호출
    fig_pca, num_proteins = run_and_plot_pca(expr_df_pca, group_map)
    
    st.plotly_chart(fig_pca, use_container_width=True)
    st.markdown(f"**Proteins used:** {num_proteins} (top {top_n} by p-adj)")

except Exception as e:
    st.error(f"PCA 시각화 중 오류가 발생했습니다: {e}")


st.markdown("---")
st.markdown("**Other visualizations:**")
col1, col2 = st.columns(2)
with col1:
    st.page_link("content/results_volcano.py", label="Volcano Plot", icon="🌋")
with col2:
    st.page_link("content/results_heatmap.py", label="Heatmap", icon="🔥")