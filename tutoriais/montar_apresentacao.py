# -*- coding: utf-8 -*-
"""Reconstroi a apresentacao NQO1 com texto e figuras EDITAVEIS (Google Slides)."""
import csv
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from PIL import Image

BASE = os.path.abspath(".")
FIG = os.path.join(BASE, "output", "figures")
TAB = os.path.join(BASE, "output", "tables")
OUT = os.path.join(BASE, "tutoriais", "APRESENTACAO_NQO1_TUTORIAIS.pptx")

TEAL = RGBColor(0x14, 0x6B, 0x63)
NAVY = RGBColor(0x1F, 0x3B, 0x5C)
ORANGE = RGBColor(0xD9, 0x5F, 0x02)
GRAY = RGBColor(0x59, 0x59, 0x59)
LIGHT = RGBColor(0xEE, 0xF4, 0xF3)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
FONT = "Arial"

W, H = Inches(13.333), Inches(7.5)
M = Inches(0.7)


def read_csv(name):
    with open(os.path.join(TAB, name), encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def br(x, dec=2):
    return f"{float(x):.{dec}f}".replace(".", ",")


# ---------------------------------------------------------------- primitivas
def txt(slide, texto, l, t, w, h, size=18, bold=False, color=NAVY, align=PP_ALIGN.LEFT,
        italic=False, spacing=1.0):
    box = slide.shapes.add_textbox(l, t, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = 0
    p = tf.paragraphs[0]
    p.alignment = align
    p.line_spacing = spacing
    r = p.add_run()
    r.text = texto
    f = r.font
    f.name, f.size, f.bold, f.italic = FONT, Pt(size), bold, italic
    f.color.rgb = color
    return box


def bullets(slide, itens, l, t, w, h, size=18, space=10):
    box = slide.shapes.add_textbox(l, t, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = 0
    for i, item in enumerate(itens):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(space)
        p.line_spacing = 1.05
        r = p.add_run()
        r.text = "\u2022  " + item
        r.font.name, r.font.size = FONT, Pt(size)
        r.font.color.rgb = NAVY


def titulo(slide, texto, sub=None):
    txt(slide, texto, M, Inches(0.42), W - 2 * M, Inches(0.75), size=30, bold=True, color=NAVY)
    bar = slide.shapes.add_shape(1, M, Inches(1.18), Inches(1.9), Pt(4))
    bar.fill.solid()
    bar.fill.fore_color.rgb = TEAL
    bar.line.fill.background()
    bar.shadow.inherit = False
    if sub:
        txt(slide, sub, M, Inches(1.32), W - 2 * M, Inches(0.4), size=14, italic=True, color=GRAY)


def rodape(slide, n, total):
    txt(slide, "NQO1 e cardiotoxicidade por doxorrubicina", M, H - Inches(0.36), Inches(7),
        Inches(0.3), size=10, color=GRAY)
    txt(slide, f"{n}/{total}", W - M - Inches(1), H - Inches(0.36), Inches(1), Inches(0.3),
        size=10, color=GRAY, align=PP_ALIGN.RIGHT)


def figura(slide, nome, l, t, w, h, legenda=None):
    caminho = os.path.join(FIG, nome)
    with Image.open(caminho) as im:
        iw, ih = im.size
    escala = min(w / iw, h / ih)
    pw, ph = int(iw * escala), int(ih * escala)
    slide.shapes.add_picture(caminho, l + int((w - pw) / 2), t + int((h - ph) / 2), pw, ph)
    if legenda:
        txt(slide, legenda, l, t + h + Inches(0.03), w, Inches(0.3), size=11,
            italic=True, color=GRAY, align=PP_ALIGN.CENTER)


def tabela(slide, cabecalho, linhas, l, t, larguras, size=13, altura=0.32):
    n_l, n_c = len(linhas) + 1, len(cabecalho)
    shp = slide.shapes.add_table(n_l, n_c, l, t, sum(larguras), Inches(altura * n_l))
    tbl = shp.table
    for i, w in enumerate(larguras):
        tbl.columns[i].width = w
    for j, h in enumerate([Inches(altura)] * n_l):
        tbl.rows[j].height = h
    for c, valor in enumerate(cabecalho):
        cell = tbl.cell(0, c)
        cell.fill.solid()
        cell.fill.fore_color.rgb = NAVY
        cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        cell.margin_left = cell.margin_right = Inches(0.06)
        p = cell.text_frame.paragraphs[0]
        r = p.add_run()
        r.text = str(valor)
        r.font.name, r.font.size, r.font.bold = FONT, Pt(size), True
        r.font.color.rgb = WHITE
    for i, linha in enumerate(linhas, start=1):
        for c, valor in enumerate(linha):
            cell = tbl.cell(i, c)
            cell.fill.solid()
            cell.fill.fore_color.rgb = LIGHT if i % 2 else WHITE
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
            cell.margin_left = cell.margin_right = Inches(0.06)
            p = cell.text_frame.paragraphs[0]
            r = p.add_run()
            r.text = str(valor)
            r.font.name, r.font.size = FONT, Pt(size)
            r.font.color.rgb = NAVY
    return tbl


def novo(prs, titulo_txt, sub=None, numero=None, total=None, notas=None):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    titulo(s, titulo_txt, sub)
    FOOTER.append(len(prs.slides))
    if notas:
        s.notes_slide.notes_text_frame.text = notas
    return s


def capa_num(prs, numero, texto, apoio):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    fundo = s.shapes.add_shape(1, 0, 0, W, H)
    fundo.fill.solid()
    fundo.fill.fore_color.rgb = NAVY
    fundo.line.fill.background()
    fundo.shadow.inherit = False
    txt(s, numero, M, Inches(2.0), Inches(2), Inches(1.6), size=80, bold=True, color=TEAL)
    txt(s, texto, M, Inches(3.3), W - 3 * M, Inches(1.2), size=40, bold=True, color=WHITE)
    txt(s, apoio, M, Inches(4.5), W - 3 * M, Inches(1.0), size=18, color=RGBColor(0xC9, 0xD8, 0xD6))
    return s


# ---------------------------------------------------------------- dados
meta = read_csv("metadados_das_amostras.csv")
grupos = {}
for row in meta:
    g = row["group"]
    grupos.setdefault(g, []).append(float(row["age"]))
total_amostras = len(meta)
n_ctrl, n_dcm, n_icm = len(grupos["CTRL"]), len(grupos["DCM"]), len(grupos["ICM"])

nqo1 = read_csv("expressao_nqo1_diferencial.csv")[0]
gene_sets = {r["conjunto"]: int(r["n"]) for r in read_csv("genes_dataset_vs_ampk.csv")}
gsea = read_csv("gsea_via_ampk.csv")[0]
gsva = read_csv("gsva_via_ampk_teste.csv")[0]
gsva_res = {r["group"]: r for r in read_csv("gsva_via_ampk_resumo_grupos.csv")}
hubs = read_csv("ppi_ampk_hubs.csv")
svac = {r["metrica"]: r for r in read_csv("SVA_comparacao.csv")}
umap = {r["group"]: r for r in read_csv("UMAP_resumo_grupos.csv")}
via = read_csv("expressao_diferencial_via_ampk.csv")
via_signif = [r for r in via if float(r["padj"]) < 0.05]
via_up = sum(1 for r in via_signif if float(r["log2FoldChange"]) > 0)
via_down = len(via_signif) - via_up

prs = Presentation()
prs.slide_width, prs.slide_height = W, H
TOTAL = 31
FOOTER = []

# 1 ----------------------------------------------------------------- capa
s = prs.slides.add_slide(prs.slide_layouts[6])
faixa = s.shapes.add_shape(1, 0, 0, W, Inches(0.32))
faixa.fill.solid()
faixa.fill.fore_color.rgb = TEAL
faixa.line.fill.background()
faixa.shadow.inherit = False
txt(s, "NQO1, doxorrubicina e cardiotoxicidade", M, Inches(1.5), W - 2 * M, Inches(0.7),
    size=40, bold=True, color=NAVY)
txt(s, "na cardiomiopatia dilatada", M, Inches(2.2), W - 2 * M, Inches(0.7), size=34,
    bold=True, color=TEAL)
linha = s.shapes.add_shape(1, M, Inches(3.05), Inches(3.2), Pt(3))
linha.fill.solid()
linha.fill.fore_color.rgb = ORANGE
linha.line.fill.background()
linha.shadow.inherit = False
txt(s, "Análise transcriptômica reprodutível em R / Bioconductor", M, Inches(3.3),
    W - 2 * M, Inches(0.5), size=20, color=GRAY)
txt(s, "L. P. Souza  ·  R. P. Santos  ·  P. A. Borges", M, Inches(5.4), W - 2 * M,
    Inches(0.4), size=16, bold=True, color=NAVY)
txt(s, "Dados públicos: GEO GSE116250  |  Via AMPK: KEGG hsa04152  |  2026", M, Inches(5.85),
    W - 2 * M, Inches(0.4), size=13, color=GRAY)
s.notes_slide.notes_text_frame.text = (
    "Apresentar o objetivo: mostrar como um pipeline reprodutível em R identifica o NQO1 "
    "downregulado na cardiomiopatia dilatada e o que isso significa para a cardiotoxicidade "
    "da doxorrubicina. Reforçar que todos os números vêm de scripts versionados.")

# 2 ----------------------------------------------------------------- roteiro
s = novo(prs, "Roteiro da aula", numero=2, total=TOTAL, notas=(
    "Cinco blocos: contexto clínico, dados e métodos, resultados gene a gene, "
    "resultados de sistema (vias e redes) e robustez/conclusões. Cada bloco "
    "corresponde a um ou dois scripts da pasta tutoriais/."))
tabela(s, ["Bloco", "Pergunta que responde", "Script"],
       [["1. Contexto", "Por que NQO1 importa na cardiotoxicidade?", "—"],
        ["2. Dados e métodos", "De onde vêm os dados e como são analisados?", "00"],
        ["3. Exploração", "As amostras se separam por grupo?", "01_UMAP"],
        ["4. Expressão diferencial", "Quais genes mudam entre DCM e controle?", "02_volcano"],
        ["5. Do gene ao sistema", "A via AMPK, a rede PPI e os escores por amostra", "03 a 06"],
        ["6. Robustez e conclusão", "O achado resiste a confundidores?", "06"]],
       M, Inches(1.85), [Inches(3.3), Inches(6.3), Inches(2.3)], size=14, altura=0.62)

# 3 ----------------------------------------------------------------- seção 1
capa_num(prs, "01", "Contexto clínico", "Doxorrubicina, estresse oxidativo e o papel do NQO1")

# 4 ----------------------------------------------------------------- DOX
s = novo(prs, "Doxorrubicina: eficácia e preço cardíaco", numero=4, total=TOTAL, notas=(
    "A antraciclina é um dos quimioterápicos mais usados e um dos mais cardiotóxicos. "
    "A toxicidade é dose-dependente e pode aparecer anos depois. "
    "Explicar que o coração tem pouca capacidade de regeneração e alta demanda energética."))
bullets(s, [
    "Antraciclina essencial no tratamento de leucemias, linfomas e tumores sólidos.",
    "Cardiotoxicidade dose-dependente: disfunção ventricular, insuficiência cardíaca e arritmias.",
    "Pode surgir durante o tratamento ou anos depois (toxicidade tardia).",
    "O miocárdio tem alta demanda energética e baixa capacidade regenerativa.",
    "Isso torna o coração especialmente sensível a dano mitocondrial e estresse oxidativo.",
], M, Inches(1.9), Inches(7.0), Inches(4.8), size=19)

# 5 ----------------------------------------------------------------- mecanismos
s = novo(prs, "Como a doxorrubicina lesiona o cardiomiócito", numero=5, total=TOTAL, notas=(
    "Quatro mecanismos que se somam. O ponto comum é a geração de espécies reativas "
    "e a disfunção mitocondrial. O NQO1 entra justamente na defesa antioxidante."))
tabela(s, ["Mecanismo", "O que acontece"],
       [["Estresse oxidativo", "Redução da quinona gera radicais livres e peroxidação lipídica"],
        ["Topoisomerase II-beta", "Inibição no cardiomiócito altera a transcrição mitocondrial"],
        ["Disfunção mitocondrial", "Perda de potencial de membrana, queda de ATP e liberação de citocromo c"],
        ["Ferroptose e apoptose", "Sobrecarga de ferro, lipoperoxidação e morte celular programada"]],
       M, Inches(1.85), [Inches(3.6), Inches(8.3)], size=14, altura=0.72)
txt(s, "Consequência comum: menos defesa antioxidante e mais vulnerabilidade ao dano.",
    M, Inches(5.6), Inches(11.5), Inches(0.5), size=16, bold=True, color=TEAL)

# 5a ---------------------------------------------------------------- Kong Fig 1
s = novo(prs, "O ciclo redox da doxorrubicina", sub="Kong et al., 2022 — Fig. 1",
         notas=("A figura tem cinco painéis. A e B: ciclo redox da quinona com NOX2/NOX4, "
                "xantina oxidase e eNOS/iNOS gera superóxido, peroxinitrito e radical hidroxila. "
                "C: no retículo sarcoplasmático, ROS ativam CaMKII. D: inibição de topoisomerase "
                "II-beta leva a quebras de mtDNA e queda de ATP. E: bloqueio do fluxo autofágico "
                "pela falha de acidificação lisossomal."))
figura(s, "artigos/Kong_2022_Fig1_estresse_oxidativo.png", M, Inches(1.45),
       Inches(4.9), Inches(5.25), "Kong et al., 2022, Int J Biol Sci")
bullets(s, [
    "A · A quinona da DOX entra em ciclo redox com NOX2, NOX4 e xantina oxidase.",
    "A · Surge a semiquinona, que reduz O2 e produz superóxido (O2\u2022-).",
    "B · A SOD converte O2\u2022- em H2O2; com Fe2+ (Fenton) vem o radical hidroxila.",
    "B · Óxido nítrico + superóxido = peroxinitrito (ONOO-).",
    "C/D · ROS ativam CaMKII no retículo; a inibição de topo-II\u03b2 quebra o mtDNA.",
    "E · A DOX bloqueia a acidificação lisossomal e trava o fluxo autofágico.",
], M + Inches(5.25), Inches(1.6), Inches(6.6), Inches(5.0), size=14, space=9)

# 5b ---------------------------------------------------------------- Kong Fig 2
s = novo(prs, "Três mortes clássicas induzidas pela DOX",
         sub="Apoptose, piroptose e ferroptose — Kong et al., 2022, Fig. 2",
         notas=("Painel A: apoptose, com PARP e p53 acionando caspases. Painel B: piroptose via "
                "TLR4/NLRP3 e caspase-1, com GSDMD/GSDME formando poros. Painel C: ferroptose, "
                "dependente de Fe2+ e peroxidação lipídica, contida por GPX4/system Xc- e GSH."))
figura(s, "artigos/Kong_2022_Fig2_morte_celular.png", M, Inches(1.45),
       Inches(11.9), Inches(3.85))
bullets(s, [
    "A · Apoptose: quebras de DNA ativam PARP e p53, que acionam caspases.",
    "B · Piroptose: TLR4 leva ao inflamassoma NLRP3 e à caspase-1; GSDMD abre poros na membrana.",
    "C · Ferroptose: Fe2+ e peroxidação lipídica; GPX4, system Xc- e GSH tentam conter.",
], M, Inches(5.45), Inches(11.9), Inches(1.7), size=14, space=7)

# 5c ---------------------------------------------------------------- Rodrigues Fig 1
s = novo(prs, "As 25 formas de morte celular regulada",
         sub="Rodrigues De Souza et al., 2025 — Fig. 1",
         notas=("Tabela com as 25 formas de morte celular regulada e a definição operacional de cada "
                "uma. Útil para mostrar que apoptose, necroptose, piroptose, ferroptose, cobreptose, "
                "oxeiptose e outras são categorias distintas, com mecanismos próprios."))
figura(s, "artigos/Rodrigues_2025_Fig1_25_mortes_celulares.png", M, Inches(1.38),
       Inches(11.9), Inches(5.1),
       "RCD = regulated cell death: forma de morte e definição operacional")

# 5d ---------------------------------------------------------------- Rodrigues Fig 3
s = novo(prs, "Multioptose: as vias conversam entre si",
         sub="Rodrigues De Souza et al., 2025 — Fig. 3",
         notas=("O modelo integrado calcula escores das 25 vias de morte por amostra e usa isso "
                "para buscar biomarcadores e alvos terapêuticos em pan-cancer. A mensagem central: "
                "as vias se comunicam e bloquear uma pode desviar a célula para outra."))
figura(s, "artigos/Rodrigues_2025_Fig3_workflow_multioptose.png", M, Inches(1.45),
       Inches(7.3), Inches(5.1), "Fluxo do modelo integrado de multioptose")
bullets(s, [
    "O modelo soma escores das 25 vias de morte, amostra por amostra.",
    "Integra camadas de multi-ômica para achar biomarcadores e alvos.",
    "As vias se comunicam: bloquear uma pode desviar a célula para outra.",
], M + Inches(7.6), Inches(1.9), Inches(4.7), Inches(3.2), size=15, space=10)

# 5e ---------------------------------------------------------------- NQO1 na historia
s = novo(prs, "Onde o NQO1 entra nessa história",
         sub="A enzima corta o ciclo redox na origem",
         notas=("Amarrar as figuras dos artigos ao achado do estudo: o NQO1 reduz a quinona em dois "
                "elétrons, evitando a semiquinona que alimenta o ciclo redox. Com NQO1 baixo na DCM, "
                "a mesma dose de DOX encontra menos defesa."))
tabela(s, ["Etapa", "Com DOX e pouco NQO1", "Com NQO1 ativo"],
       [["Ciclo redox", "Semiquinona recicla e gera superóxido sem parar",
         "Redução em dois elétrons: sem semiquinona"],
        ["Radicais", "H2O2, \u2022OH e ONOO- atacam proteína, lipídio e DNA",
         "Menos ROS e menos peroxidação lipídica"],
        ["Ferroptose", "Fe2+ + peroxidação lipídica levam à morte por ferro",
         "Menos peroxidação: menos sinal ferroptótico"],
        ["Resultado", "Morte do cardiomiócito", "Mais resistência ao dano"]],
       M, Inches(1.8), [Inches(2.5), Inches(4.9), Inches(4.5)], size=13, altura=0.66)
txt(s, "Na DCM o NQO1 cai: a defesa antioxidante diminui justamente onde a DOX age.",
    M, Inches(5.3), Inches(11.5), Inches(0.9), size=17, bold=True, color=TEAL)

# 6 ----------------------------------------------------------------- NQO1
s = novo(prs, "NQO1: a enzima que detoxifica quinonas", numero=6, total=TOTAL, notas=(
    "NQO1 é NAD(P)H quinona oxidorredutase 1. Faz redução de dois elétrons, evitando a "
    "formação de semiquinonas reativas. É induzida por Nrf2 via elemento de resposta "
    "antioxidante. Na prática: é uma enzima citoprotetora."))
bullets(s, [
    "NAD(P)H quinona oxidorredutase 1: reduz quinonas por transferência de dois elétrons.",
    "Evita a formação de semiquinonas reativas e de espécies reativas de oxigênio.",
    "Induzida por Nrf2 pelo elemento de resposta antioxidante (ARE).",
    "Citoprotetora: protege o cardiomiócito do dano oxidativo induzido pela DOX.",
    "Hipótese do estudo: na cardiomiopatia dilatada, o NQO1 está downregulado.",
], M, Inches(1.9), Inches(7.0), Inches(4.8), size=19)

# 7 ----------------------------------------------------------------- pergunta
s = novo(prs, "Pergunta e desenho do estudo", numero=7, total=TOTAL, notas=(
    "Pergunta principal e três perguntas secundárias. O desenho é observacional, "
    "com dados públicos de transcriptoma de ventrículo esquerdo. Deixar claro que "
    "não há intervenção nem seguimento clínico."))
txt(s, "O NQO1 está diferencialmente expresso na cardiomiopatia dilatada e a via "
       "AMPK está perturbada nesse contexto?",
    M, Inches(1.8), Inches(11.9), Inches(1.0), size=22, bold=True, color=TEAL)
tabela(s, ["Objetivo específico", "Abordagem"],
       [["Quantificar o NQO1 entre DCM e controle", "limma, ajustado por sexo e idade"],
        ["Testar a via AMPK como conjunto", "GSEA e GSVA (KEGG hsa04152)"],
        ["Ver o sistema, não só o gene", "Rede PPI no STRING (score >= 700)"],
        ["Checar confundidores técnicos", "SVA para fatores latentes e xCell para composição celular"]],
       M, Inches(3.05), [Inches(5.6), Inches(6.3)], size=14, altura=0.6)

# 8 ----------------------------------------------------------------- seção 2
capa_num(prs, "02", "Dados e métodos", f"{total_amostras} amostras de ventrículo esquerdo, pipeline reprodutível")

# 9 ----------------------------------------------------------------- dataset
s = novo(prs, "Dataset GSE116250", numero=9, total=TOTAL, notas=(
    "RNA-seq de ventrículo esquerdo humano. Usamos só DCM e controle no teste "
    "estatístico; o grupo ICM entra apenas como contexto no UMAP. "
    "Idade e sexo entram como covariáveis."))
tabela(s, ["Grupo", "n", "Idade média (anos)", "Papel na análise"],
       [["CTRL", str(n_ctrl), br(sum(grupos["CTRL"]) / n_ctrl, 1), "Controle usado no teste"],
        ["DCM", str(n_dcm), br(sum(grupos["DCM"]) / n_dcm, 1), "Caso usado no teste"],
        ["ICM", str(n_icm), br(sum(grupos["ICM"]) / n_icm, 1), "Só contexto no UMAP"],
        ["Total", str(total_amostras), br(sum(sum(v) for v in grupos.values()) / total_amostras, 1), "—"]],
       M, Inches(1.85), [Inches(2.2), Inches(1.6), Inches(3.2), Inches(4.9)], size=15, altura=0.6)
txt(s, "Expressão em RPKM · 18.526 genes após filtro de expressão",
    M, Inches(4.9), Inches(11.5), Inches(0.4), size=15, bold=True, color=TEAL)

# 10 ---------------------------------------------------------------- pipeline
s = novo(prs, "Pipeline: sete etapas, sete scripts", numero=10, total=TOTAL, notas=(
    "A pasta tutoriais/ tem um script por etapa. O script 00 prepara os dados uma única "
    "vez e guarda cache; os demais reaproveitam. A ordem dos arquivos é a ordem de execução."))
tabela(s, ["Script", "O que faz", "Saída principal"],
       [["00", "Configuração, download do GEO e da via KEGG", "cache em data/processed"],
        ["01_UMAP", "Exploração não supervisionada das amostras", "UMAP_amostras.png"],
        ["02_volcano", "limma DCM vs CTRL + volcano", "expressao_diferencial_completa.csv"],
        ["03_gene_sets", "Genes do dataset x genes da via AMPK", "Venn_dataset_AMPK.png"],
        ["04_rede_PPI", "Rede de interações no STRING", "rede_ppi_ampk.png"],
        ["05_extras_GSEA_GSVA", "GSEA e GSVA da via AMPK", "GSEA e GSVA (png/csv)"],
        ["06_batch_deconvolution", "SVA e xCell", "SVA_comparacao.csv, xCell_teste.csv"]],
       M, Inches(1.8), [Inches(2.6), Inches(5.9), Inches(3.9)], size=13, altura=0.52)
txt(s, "Tudo reprodutível: Rscript executar_tutoriais.R",
    M, Inches(6.0), Inches(11.5), Inches(0.4), size=14, bold=True, color=TEAL)

# 11 ---------------------------------------------------------------- limma
s = novo(prs, "Modelo estatístico (limma)", numero=11, total=TOTAL, notas=(
    "limma ajusta um modelo linear por gene com efeitos de grupo, sexo e idade. "
    "A moderação empírica da variância estabiliza as estimativas com poucas amostras. "
    "O ajuste de múltiplas comparações é Benjamini-Hochberg."))
bullets(s, [
    "Modelo linear por gene: grupo (DCM vs CTRL) + sexo + idade.",
    "Moderação empírica da variância: estabiliza o teste com amostras pequenas.",
    "log2FC negativo indica menor expressão na DCM.",
    "Controle de falsos positivos por Benjamini-Hochberg (padj).",
    "Corte de significância: padj < 0,05; efeito relevante: |log2FC| >= 0,585.",
], M, Inches(1.9), Inches(7.2), Inches(4.6), size=18)

# 12 ---------------------------------------------------------------- seção 3
capa_num(prs, "03", "Resultados: exploração", "As amostras se organizam por grupo antes de qualquer teste?")

# 13 ---------------------------------------------------------------- UMAP
s = novo(prs, "UMAP: os controles formam um grupo à parte", numero=13, total=TOTAL,
         notas=("UMAP não supervisionado, colorido depois pelos rótulos. CTRL separa; DCM e ICM "
                "se sobrepõem, o que é esperado: são doenças diferentes com fenótipo final parecido. "
                "Eixos do UMAP não têm interpretação biológica direta."))
figura(s, "UMAP_amostras.png", M, Inches(1.5), Inches(7.6), Inches(5.3),
       "Projeção UMAP das amostras (n = 68)")

# 14 ---------------------------------------------------------------- UMAP tabela
s = novo(prs, "UMAP em números", numero=14, total=TOTAL,
         notas="Médias e desvios das coordenadas por grupo. O CTRL tem menor dispersão e centro distante.")
tabela(s, ["Grupo", "n", "UMAP1 (média)", "UMAP2 (média)", "DP UMAP1"],
       [[g, umap[g]["n"], br(umap[g]["UMAP1_medio"]), br(umap[g]["UMAP2_medio"]),
         br(umap[g]["sd_UMAP1"])] for g in ("CTRL", "DCM", "ICM")],
       M, Inches(1.9), [Inches(2.0), Inches(1.2), Inches(2.6), Inches(2.6), Inches(2.4)],
       size=15, altura=0.55)
txt(s, "O resultado depende de n_neighbors, min_dist e do número de genes usados.",
    M, Inches(4.6), Inches(11.5), Inches(0.4), size=14, italic=True, color=GRAY)

# 15 ---------------------------------------------------------------- volcano/panorama
s = novo(prs, "Expressão diferencial: o panorama", numero=15, total=TOTAL,
         notas=("Quase metade dos genes testados tem padj < 0,05, o que mostra que a DCM é um "
                "estado transcriptômico muito alterado — e por isso o efeito precisa ser lido "
                "junto com o tamanho de efeito."))
tabela(s, ["Métrica", "Valor"],
       [["Genes testados", "18.526"],
        ["DEGs com padj < 0,05", svac["DEGs padj < 0,05"]["sem_SVA"]],
        ["DEGs com padj < 0,05 e |log2FC| >= 0,585", svac["DEGs padj < 0,05 e |log2FC| >= 0,585"]["sem_SVA"]],
        ["Genes da via AMPK testados", str(len(via))],
        ["Genes da via AMPK com padj < 0,05", f"{len(via_signif)} ({via_up} up, {via_down} down)"]],
       M, Inches(1.85), [Inches(6.6), Inches(4.6)], size=15, altura=0.6)
txt(s, "A DCM é um estado transcriptômico muito alterado: o tamanho de efeito importa.",
    M, Inches(5.9), Inches(11.5), Inches(0.4), size=14, italic=True, color=GRAY)

# 16 ---------------------------------------------------------------- NQO1 resultado
s = novo(prs, "NQO1 está downregulado na DCM", numero=16, total=TOTAL,
         notas=("Resultado central do estudo. log2FC -1,66 corresponde a cerca de um terço da "
                "expressão do controle. padj 6,5e-13 é muito robusto."))
tabela(s, ["Gene", "log2FC", "Expressão média", "p-valor", "padj", "Direção"],
       [[nqo1["gene"], br(nqo1["log2FoldChange"]), br(nqo1["AveExpr"]),
         "1,4e-16", "6,5e-13", "down na DCM"]],
       M, Inches(1.9), [Inches(1.7), Inches(1.9), Inches(2.6), Inches(1.9), Inches(1.9), Inches(2.2)],
       size=15, altura=0.55)
txt(s, "-1,66 em log2 = aproximadamente um terço da expressão do controle.",
    M, Inches(3.2), Inches(11.5), Inches(0.45), size=18, bold=True, color=TEAL)
bullets(s, [
    "Menos NQO1 significa menos detoxificação de quinonas no miocárdio.",
    "Um coração já doente fica mais vulnerável ao dano da doxorrubicina.",
    "Isso conecta o achado transcriptômico à hipótese clínica de susceptibilidade.",
], M, Inches(3.95), Inches(11.6), Inches(2.0), size=18)

# 17 ---------------------------------------------------------------- seção 4
capa_num(prs, "04", "Do gene ao sistema", "Gene sets, via AMPK e rede de interações proteicas")

# 18 ---------------------------------------------------------------- gene sets
s = novo(prs, "Quais genes entraram na análise", numero=18, total=TOTAL,
         notas=("Distinção importante: estar no dataset não é o mesmo que pertencer à via. "
                "A interseção é o que o limma, a rede e o GSEA efetivamente usam."))
figura(s, "Venn_dataset_AMPK.png", M, Inches(1.5), Inches(6.6), Inches(4.7),
       "Dataset filtrado x via AMPK (KEGG hsa04152)")
tabela(s, ["Conjunto", "n"],
       [[k, v] for k, v in gene_sets.items()],
       M + Inches(7.0), Inches(1.9), [Inches(3.4), Inches(1.6)], size=13, altura=0.55)

# 19 ---------------------------------------------------------------- rede PPI
s = novo(prs, "Rede PPI da via AMPK no STRING", numero=19, total=TOTAL,
         notas=("Rede com 105 proteínas e 837 interações de alta confiança. Tamanho do nó = grau; "
                "cor = log2FC. Hubs concentram a comunicação da via."))
figura(s, "rede_ppi_ampk.png", M, Inches(1.45), Inches(8.4), Inches(5.1),
       "Rede PPI (score >= 700) · 105 nós e 837 arestas")
bullets(s, [
    "Grafo: nós = proteínas, arestas = interações.",
    "Score >= 700 de 1.000 = alta confiança.",
    "Hubs = proteínas mais conectadas.",
], M + Inches(8.7), Inches(2.0), Inches(3.6), Inches(2.4), size=15)

# 20 ---------------------------------------------------------------- hubs
s = novo(prs, "Os hubs da rede", numero=20, total=TOTAL,
         notas=("AKT1 e MTOR lideram. São nós de sinalização de crescimento e metabolismo, coerentes "
                "com o papel energético da AMPK no coração."))
tabela(s, ["Gene", "Grau", "log2FC"],
       [[h["gene"], h["degree"], br(h["log2FoldChange"])] for h in hubs[:10]],
       M, Inches(1.8), [Inches(2.6), Inches(2.2), Inches(2.6)], size=14, altura=0.44)
txt(s, "Grau = número de conexões na rede. AKT1 e MTOR concentram a sinalização.",
    M + Inches(6.2), Inches(2.1), Inches(5.8), Inches(1.2), size=16, color=GRAY)

# 21 ---------------------------------------------------------------- seção 5
capa_num(prs, "05", "Vias e robustez", "GSEA, GSVA, fatores latentes e composição celular")

# 22 ---------------------------------------------------------------- GSEA
s = novo(prs, "GSEA: a via AMPK não enriquece", numero=22, total=TOTAL,
         notas=("NES negativo e p 0,31: não há enriquecimento significativo do conjunto neste recorte. "
                "Resultado negativo também é resultado: a via não se comporta como bloco."))
tabela(s, ["Métrica", "Valor"],
       [["Genes da via no ranking", gsea["n_genes"]],
        ["NES", br(gsea["NES"])],
        ["p-valor", br(gsea["p_value"], 3)]],
       W - M - Inches(4.4), Inches(1.8), [Inches(2.4), Inches(2.0)], size=14, altura=0.55)
txt(s, "NES = -1,06 e p = 0,31: sem enriquecimento significativo.",
    M, Inches(1.85), Inches(7.4), Inches(0.5), size=18, bold=True, color=TEAL)
figura(s, "GSEA_AMPK_enrichment_plot.png", M, Inches(2.5), Inches(11.9), Inches(4.3),
       "Curva de enriquecimento da via AMPK")
txt(s, "A via como bloco não se move — mas genes individuais, sim.",
    M, Inches(6.85), Inches(11.5), Inches(0.35), size=12, italic=True, color=GRAY)

# 23 ---------------------------------------------------------------- GSVA
s = novo(prs, "GSVA: escore da via por amostra", numero=23, total=TOTAL,
         notas=("GSVA transforma a via em um escore por amostra. Aqui a diferença é significativa "
                "por Wilcoxon, com escores maiores no controle."))
tabela(s, ["Grupo", "n", "Média", "Mediana"],
       [[g, gsva_res[g]["n"], br(gsva_res[g]["media"], 3), br(gsva_res[g]["mediana"], 3)]
        for g in ("CTRL", "DCM")],
       W - M - Inches(4.6), Inches(1.8), [Inches(1.4), Inches(1.0), Inches(1.3), Inches(1.6)],
       size=13, altura=0.5)
txt(s, f"Wilcoxon p = {br(gsva['p_value'], 5)}  ·  -log10(p) = {br(gsva['neg_log10_p'])}",
    M, Inches(1.85), Inches(7.2), Inches(0.5), size=18, bold=True, color=TEAL)
figura(s, "GSVA_AMPK_boxplot.png", M, Inches(2.5), Inches(11.9), Inches(4.2),
       "Escore GSVA da via AMPK por grupo")
txt(s, "O escore médio cai de 0,115 (CTRL) para -0,058 (DCM).",
    M, Inches(6.85), Inches(11.5), Inches(0.35), size=12, italic=True, color=GRAY)

# 24 ---------------------------------------------------------------- SVA
s = novo(prs, "Robustez: correção de fatores latentes (SVA)", numero=24, total=TOTAL,
         notas=("SVA estima fatores latentes sem usar os rótulos. Com dois fatores no modelo, o efeito "
                "do NQO1 praticamente não muda. Isso sustenta o achado contra confundidores técnicos."))
tabela(s, ["Métrica", "Sem SVA", "Com SVA"],
       [["NQO1 log2FC", br(svac["NQO1 log2FC"]["sem_SVA"]), br(svac["NQO1 log2FC"]["com_SVA"])],
        ["NQO1 padj", "6,5e-13", "7,4e-13"],
        ["DEGs padj < 0,05", svac["DEGs padj < 0,05"]["sem_SVA"], svac["DEGs padj < 0,05"]["com_SVA"]],
        ["Correlação de log2FC", "1,000", br(svac["Correlação log2FC (todos os genes)"]["com_SVA"], 3)]],
       M, Inches(1.8), [Inches(3.6), Inches(1.9), Inches(1.9)], size=14, altura=0.55)
put = s.shapes.add_picture(os.path.join(FIG, "SVA_logFC_antes_depois.png"),
                           M + Inches(8.0), Inches(1.75), Inches(3.9))
txt(s, "O efeito do NQO1 sobrevive à correção: log2FC -1,66 -> -1,65.",
    M, Inches(5.0), Inches(7.6), Inches(0.9), size=17, bold=True, color=TEAL)

# 25 ---------------------------------------------------------------- xCell heatmap
s = novo(prs, "Deconvolução celular (xCell)", numero=25, total=TOTAL,
         notas=("Bulk de tecido é mistura de tipos celulares. xCell estima escores relativos. "
                "Serve para checar se a diferença de NQO1 é só composição celular."))
figura(s, "xCell_heatmap.png", M, Inches(1.5), Inches(6.2), Inches(4.9),
       "Escores xCell por amostra")
bullets(s, [
    "Bulk cardíaco = mistura de tipos celulares.",
    "xCell estima escores relativos, não frações absolutas.",
    "Ajuda a descartar que a diferença seja só composição.",
], M + Inches(6.6), Inches(1.9), Inches(5.9), Inches(2.6), size=16)

# 26 ---------------------------------------------------------------- xCell boxplots
s = novo(prs, "xCell: tipos celulares que mudam", numero=26, total=TOTAL,
         notas=("Miocitos e endotelio linfatico sobem; macrofagos M1 e M2 caem. Ler como escores "
                "relativos, com cautela, porque a composição influencia todos os escores."))
figura(s, "xCell_boxplots.png", M, Inches(1.5), Inches(8.0), Inches(5.2),
       "Escores por tipo celular, DCM vs CTRL")
bullets(s, [
    "Miocitos e endotelio linfatico aumentados.",
    "Macrofagos M1 e M2 reduzidos.",
    "Interpretar como escores relativos.",
], M + Inches(8.3), Inches(2.0), Inches(4.0), Inches(2.6), size=15)

# 27 ---------------------------------------------------------------- integração
s = novo(prs, "O que os resultados contam juntos", numero=27, total=TOTAL,
         notas=("Fechar a narrativa: NQO1 baixo, via AMPK sem enriquecimento global mas com "
                "heterogeneidade, hubs de sinalização preservados e composição celular alterada."))
bullets(s, [
    "NQO1 está reduzido de forma robusta na DCM (log2FC -1,66; padj 6,5e-13).",
    "A via AMPK não enriquece como bloco no GSEA, mas varia no escore por amostra (GSVA).",
    "A rede PPI mantém hubs de sinalização (AKT1, MTOR, PRKAA1/2).",
    "A composição celular muda, e o achado do NQO1 resiste à correção por SVA.",
], M, Inches(1.9), Inches(11.8), Inches(3.2), size=19)
txt(s, "Implicação: menor defesa antioxidante pode amplificar a cardiotoxicidade da DOX.",
    M, Inches(5.6), Inches(11.8), Inches(0.6), size=18, bold=True, color=ORANGE)

# 28 ---------------------------------------------------------------- limitações
s = novo(prs, "Limitações", numero=28, total=TOTAL,
         notas=("Ser honesto: dados observacionais, sem seguimento, sem dose de DOX, "
                "e bulk sem resolução de célula única."))
bullets(s, [
    "Estudo observacional com dados públicos: não há seguimento clínico.",
    "Não há informação de dose ou tempo de exposição à doxorrubicina.",
    "Transcriptoma bulk: não resolve heterogeneidade célula a célula.",
    "Validação experimental (qPCR, western blot, modelos celulares) ainda é necessária.",
], M, Inches(1.9), Inches(11.7), Inches(3.4), size=19)

# 29 ---------------------------------------------------------------- conclusões
s = novo(prs, "Conclusões", numero=29, total=TOTAL,
         notas=("Três mensagens finais: NQO1 é biomarcador candidato, a via AMPK é heterogênea "
                "e o pipeline é reprodutível e reutilizável."))
bullets(s, [
    "NQO1 downregulado na DCM é o achado central, robusto a ajuste e a fatores latentes.",
    "A via AMPK não se comporta como bloco único: genes individuais e escores importam.",
    "O pipeline em R/Bioconductor é reprodutível e pode ser reaplicado a outras coortes.",
], M, Inches(1.9), Inches(11.8), Inches(3.0), size=19)

# 30 ---------------------------------------------------------------- reprodutibilidade
s = novo(prs, "Reprodutibilidade", numero=30, total=TOTAL,
         notas=("Mostrar a pasta tutoriais e o comando único. Cada script tem comentarios "
                "didaticos e imprime FIM DO SCRIPT ao terminar."))
tabela(s, ["Arquivo", "Função"],
       [["executar_tutoriais.R", "Roda a aula inteira na ordem correta"],
        ["tutoriais/00 a 06", "Um script por etapa, com comentários didáticos"],
        ["output/tables e output/figures", "Tabelas e figuras geradas"],
        ["data/processed", "Cache dos dados e do resultado do limma"]],
       M, Inches(1.8), [Inches(5.0), Inches(6.5)], size=14, altura=0.55)
txt(s, "Rscript executar_tutoriais.R", M, Inches(4.6), Inches(11.5), Inches(0.5),
    size=20, bold=True, color=TEAL)
txt(s, "Os números desta apresentação vêm das tabelas em output/tables.",
    M, Inches(5.2), Inches(11.5), Inches(0.4), size=13, italic=True, color=GRAY)

# 31 ---------------------------------------------------------------- referências
s = novo(prs, "Referências", numero=31, total=TOTAL,
         notas="Citar as fontes dos dados e dos métodos. Ver README.md para a lista completa.")
bullets(s, [
    "GEO: GSE116250 (transcriptoma de ventrículo esquerdo humano).",
    "Kong et al. (2022) Underlying the mechanisms of doxorubicin-induced acute cardiotoxicity. Int J Biol Sci.",
    "Rodrigues De Souza et al. (2025) Integrated multi-optosis model for pan-cancer discovery.",
    "Ritchie et al. (2015) limma. Nucleic Acids Research.",
    "Szklarczyk et al. (2021) STRING. Nucleic Acids Research.",
    "Hänzelmann et al. (2013) GSVA. BMC Bioinformatics; Korotkevich et al. (2021) fgsea. bioRxiv.",
    "Leek et al. (2012) sva. Bioinformatics; Aran et al. (2017) xCell. Genome Biology.",
], M, Inches(1.8), Inches(11.8), Inches(4.4), size=14, space=7)
txt(s, "Textos completos em tutoriais/artigos/. Lista detalhada em README.md.",
    M, Inches(6.35), Inches(11.5), Inches(0.35), size=12, italic=True, color=GRAY)

# numeracao automatica das paginas: so slides de conteudo recebem rodape
total_slides = len(prs.slides)
for indice in FOOTER:
    rodape(prs.slides[indice - 1], indice, total_slides)

prs.save(OUT)
print("slides:", len(prs.slides))
print("salvo:", OUT)
