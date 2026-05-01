import streamlit as st
import time
from trp53_core import *
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.pagesizes import letter
from io import BytesIO
import matplotlib.pyplot as plt

st.set_page_config(page_title="p53 Analyzer", layout="centered")
st.title("🧬 p53 Analyzer")

HUMAN_MUT, MOUSE_MUT = ["R175H"], ["R172H"]
HUMAN_TO_MOUSE = {"R175H": "R172H"}
MOUSE_TO_HUMAN = {v: k for k, v in HUMAN_TO_MOUSE.items()}

species = st.radio("Espécie:", ["Humano", "Rato"])
mutation = st.text_input("Mutação (ex: R175H)", "R175H").upper()

if mutation:
    try:
        if species == "Humano" and mutation in MOUSE_MUT:
            st.error("Mutação de rato detetada em contexto humano.")
            st.stop()
        if species == "Rato" and mutation in HUMAN_MUT:
            st.error("Mutação humana detetada em contexto de rato.")
            st.stop()

        pos = int(mutation[1:-1])
        human_pos = map_position(pos, species == "Rato")
        equiv = HUMAN_TO_MOUSE.get(mutation) if species == "Humano" else MOUSE_TO_HUMAN.get(mutation)

        st.subheader("Mapeamento e Cromossomas")
        c_map1, c_map2 = st.columns(2)
        with c_map1:
            st.write("**Humano (TP53)**")
            st.write(f"Posição: {human_pos}")
            st.write(f"Cromossoma: {CHROMOSOME_MAP['Humano']}")
        with c_map2:
            st.write("**Rato (Trp53)**")
            st.write(f"Posição: {pos if species == 'Rato' else map_position(pos, False)}")
            st.write(f"Cromossoma: {CHROMOSOME_MAP['Rato']}")

        clinvar = get_clinvar(mutation[0], human_pos, mutation[-1])

        with st.spinner("A calcular score biológico..."):
            time.sleep(0.5)
            score, details = compute_score(mutation, human_pos, clinvar, equiv is not None)

        st.divider()
        col1, col2 = st.columns(2)
        col1.metric("Classificação", classify(score))

        with col2:
            with st.popover("ℹ️ Escala de Classificação"):
                st.write("0: Benigno")
                st.write("1 - 2: Indeterminado")
                st.write("3 - 5: Potencialmente maligno")
                st.write("≥ 6: Maligno")
            st.markdown(f"<h2 style='text-align: center;'>{score}/11</h2>", unsafe_allow_html=True)

        st.pyplot(plot_protein_map(human_pos, mutation))

        if mutation in TREATMENTS_DB:
            with st.expander("Ver tratamentos potenciais"):
                for t in TREATMENTS_DB[mutation]:
                    st.write(f"• {t}")
        else:
            st.info("Tratamentos específicos não listados para esta variante.")

        # ✅ PDF COM FORMATAÇÃO MELHORADA
        report = generate_report(mutation, pos, clinvar, score, details, species, equiv, human_pos)

        buffer = BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=letter)
        styles = getSampleStyleSheet()

        title_style = ParagraphStyle('title', parent=styles['Normal'], fontSize=16, leading=20, spaceAfter=12, fontName="Helvetica-Bold")
        section_style = ParagraphStyle('section', parent=styles['Normal'], fontSize=13, leading=16, spaceAfter=10, spaceBefore=10, fontName="Helvetica-Bold")
        normal_style = ParagraphStyle('normal', parent=styles['Normal'], fontSize=10, leading=14, spaceAfter=6)

        content = []

        for line in report.split("\n"):
            if "RELATÓRIO DE ANÁLISE" in line:
                content.append(Paragraph(line, title_style))
            elif line.startswith("---"):
                content.append(Spacer(1, 10))
                content.append(Paragraph(line.replace("---", ""), section_style))
            elif line.strip() == "":
                content.append(Spacer(1, 8))
            else:
                content.append(Paragraph(line, normal_style))

        img_buf = BytesIO()
        fig_pdf = plot_protein_map(human_pos, mutation)
        fig_pdf.savefig(img_buf, format='png', bbox_inches='tight')
        img_buf.seek(0)

        content.append(Spacer(1, 12))
        content.append(Image(img_buf, width=450, height=120))
        plt.close(fig_pdf)

        doc.build(content)
        buffer.seek(0)

        st.download_button("📥 Descarregar Relatório PDF Completo", buffer.getvalue(), f"Relatorio_p53_{mutation}.pdf")

    except ValueError:
        st.error("Formato de mutação inválido.")
    except Exception as e:
        st.error(f"Erro: {e}")