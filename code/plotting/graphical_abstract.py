"""Editable graphical abstract, with coordinates on the 1983 × 793 canvas."""

from contextlib import contextmanager
from pathlib import Path
import subprocess
import math
import xml.etree.ElementTree as ET
from PIL import ImageFont


OUTPUT = Path(__file__).resolve().parents[2] / "figures" / "graphical_abstract.svg"
WIDTH, HEIGHT = 1983, 793
SVG_NS = "http://www.w3.org/2000/svg"
INKSCAPE_NS = "http://www.inkscape.org/namespaces/inkscape"
ET.register_namespace("", SVG_NS)
ET.register_namespace("inkscape", INKSCAPE_NS)
INK = "#090765"
EDGE = "#687fb9"
GREEN = "#00c85b"


def tag(name):
    return f"{{{SVG_NS}}}{name}"


class Canvas:
    def __init__(self):
        self.root = ET.Element(tag("svg"), {
            "width": "1983", "height": "793", "viewBox": "0 0 1983 793",
            "version": "1.1", "font-family": "Arial, Helvetica, sans-serif",
        })
        ET.SubElement(self.root, tag("title")).text = "Single-cell multi-omics regulatory network inference"
        ET.SubElement(self.root, tag("desc")).text = (
            "Independent TF–gene, TF–region and region–gene fits. Exact genomic-region "
            "matching produces region-resolved networks for fetal-cortex examples, "
            "with direct TF–gene overlap assessed separately; "
            "direct TF–gene networks describe erythroid target turnover."
        )
        self.parent = self.root
        self.counter = 0
        self.defs = ET.SubElement(self.root, tag("defs"))

    def element(self, name, **attrs):
        self.counter += 1
        attrs.setdefault("id", f"object-{self.counter:03d}")
        attrs = {key.replace("_", "-"): str(value) for key, value in attrs.items() if value is not None}
        return ET.SubElement(self.parent, tag(name), attrs)

    @contextmanager
    def group(self, name, label=None, layer=False):
        previous = self.parent
        attrs = {"id": name}
        if label:
            attrs[f"{{{INKSCAPE_NS}}}label"] = label
        if layer:
            attrs[f"{{{INKSCAPE_NS}}}groupmode"] = "layer"
        self.parent = ET.SubElement(previous, tag("g"), attrs)
        yield
        self.parent = previous

    def gradient(self, name, stops, radial=False, **attrs):
        element = ET.SubElement(self.defs, tag("radialGradient" if radial else "linearGradient"),
                                {"id": name, **{k: str(v) for k, v in attrs.items()}})
        for offset, color in stops:
            ET.SubElement(element, tag("stop"), {"offset": offset, "stop-color": color})

    def rect(self, x, y, width, height, fill, radius=0, stroke=None, sw=1, **attrs):
        return self.element("rect", x=x, y=y, width=width, height=height, rx=radius,
                            fill=fill, stroke=stroke, stroke_width=sw if stroke else None, **attrs)

    def circle(self, x, y, radius, fill, stroke=None, sw=1, **attrs):
        return self.element("circle", cx=x, cy=y, r=radius, fill=fill,
                            stroke=stroke, stroke_width=sw if stroke else None, **attrs)

    def path(self, d, color=EDGE, width=3.4, fill="none", **attrs):
        return self.element("path", d=d, fill=fill, stroke=color,
                            stroke_width=width if color else None,
                            stroke_linejoin="round", stroke_linecap="round", **attrs)

    def polygon(self, points, fill, stroke=None, sw=1, **attrs):
        return self.element("polygon", points=" ".join(f"{x:g},{y:g}" for x, y in points),
                            fill=fill, stroke=stroke, stroke_width=sw if stroke else None, **attrs)

    def text(self, x, y, value, size=28, weight="normal", anchor="middle", color=INK,
             italic=False, length=None, **attrs):
        transform = None
        if length is not None:
            font_name = "Arial Bold.ttf" if weight == "bold" else "Arial Italic.ttf" if italic else "Arial.ttf"
            font = ImageFont.truetype(font_name, 1000)
            natural_width = font.getlength(value) * size / 1000
            transform = f"translate({x:g},{y:g}) scale({length / natural_width:.8f},1)"
            x, y = 0, 0
        element = self.element("text", x=x, y=y, fill=color, font_size=size,
                               font_weight=weight, text_anchor=anchor,
                               font_style="italic" if italic else None,
                               transform=transform,
                               **attrs)
        element.text = value
        return element

    def multiline(self, x, y, lines, size=28, spacing=29, **kwargs):
        element = self.text(x, y, "", size=size, **kwargs)
        for index, line in enumerate(lines):
            ET.SubElement(element, tag("tspan"), {"x": str(x), "y": str(y + spacing * index)}).text = line
        return element

    def arrow(self, x1, y1, x2, y2, color=EDGE, width=3.4, head=18, dash=None):
        dx, dy = x2 - x1, y2 - y1
        distance = math.hypot(dx, dy)
        ux, uy = dx / distance, dy / distance
        bx, by = x2 - head * ux, y2 - head * uy
        self.path(f"M{x1:g},{y1:g} L{bx:g},{by:g}", color, width, stroke_dasharray=dash)
        self.polygon([(x2, y2), (bx - uy * head * .5, by + ux * head * .5),
                      (bx + uy * head * .5, by - ux * head * .5)], color)

    def diamond(self, x, y, rx=24, ry=24):
        self.polygon([(x, y - ry), (x + rx, y), (x, y + ry), (x - rx, y)],
                     "url(#region-node)", "#f8ffff", .8)

    def tf(self, x, y, width=28, height=None):
        self.rect(x - width / 2, y - (height or width) / 2, width, height or width,
                  "url(#tf-node)", stroke="#eef6fb", sw=.8)

    def gene(self, x, y, radius=22):
        self.circle(x, y, radius, "url(#gene-node)", "#fff9ed", 1.2)


RNA_COLORS = [
    ["#b0d5fd", "#a1d1fd", "#d8e7fd", "#61b4fd", "#019afd", "#86cbfe", "#019afc"],
    ["#f5f9fd", "#85bffd", "#0280fd", "#f6fafe", "#85c3fe", "#cde7fd", "#88ccfc"],
    ["#dcebfd", "#fcfdfe", "#90c4fd", "#e5f0fd", "#018bfd", "#0096fc", "#89cdfd"],
    ["#bfdcfd", "#88c0fd", "#dce3fd", "#6eb9fd", "#deedfe", "#60c0fe", "#91d3fe"],
]
ATAC_COLORS = [
    ["#d9efd8", "#d6eed6", "#8dd295", "#d3ecd3", "#63c473", "#7ccb93", "#82cd9b"],
    ["#e4f4e2", "#a8dcb2", "#ddf2de", "#6ec87e", "#d8eed7", "#018354", "#97d4af"],
    ["#afddb4", "#e2f2e1", "#bde6c1", "#35be5d", "#69c480", "#d5edd7", "#9bd8ac"],
    ["#ddf1db", "#aadbac", "#dcf0db", "#bfe4c2", "#02b540", "#9fdbb1", "#02ba3e"],
]


def definitions(c):
    c.gradient("input-panel", [("0%", "#fafafb"), ("50%", "#f4f4f7"), ("100%", "#e9eaf0")],
               radial=True, cx="50%", cy="46%", r="90%")
    c.gradient("method-panel", [("0%", "#e9eaf0"), ("52%", "#eff0f4"), ("100%", "#e9eaf0")],
               x1="0%", y1="0%", x2="100%", y2="100%")
    c.gradient("application-panel", [("0%", "#eeeff3"), ("100%", "#e9eaf0")],
               x1="0%", y1="0%", x2="100%", y2="100%")
    c.gradient("caption-fill", [("0%", "#e2e3ea"), ("100%", "#dcdee7")],
               x1="0%", y1="0%", x2="100%", y2="100%")
    c.gradient("tf-node", [("0%", "#2d2d2b"), ("100%", "#171817")],
               x1="0%", y1="0%", x2="100%", y2="100%")
    c.gradient("region-node", [("0%", "#07c76f"), ("100%", "#00b94d")],
               x1="0%", y1="0%", x2="20%", y2="100%")
    c.gradient("gene-node", [("0%", "#ff8434"), ("100%", "#ff6c13")],
               radial=True, cx="35%", cy="35%", r="75%")
    c.gradient("cell-body", [("0%", "#fff0f6"), ("55%", "#ffc0d8"), ("100%", "#ff8fb8")],
               radial=True, cx="37%", cy="30%", r="76%")
    c.gradient("cell-nucleus", [("0%", "#8248a3"), ("100%", "#54216f")],
               radial=True, cx="33%", cy="32%", r="76%")
    c.gradient("neuron-blue", [("0%", "#79dfff"), ("100%", "#00a8f0")],
               radial=True, cx="36%", cy="30%", r="73%")
    c.gradient("neuron-green", [("0%", "#95e1ad"), ("100%", "#3cbd75")],
               radial=True, cx="36%", cy="30%", r="73%")


def matrix(c, name, x, colors):
    with c.group(name, f"{name.upper()} matrix cells"):
        for row, values in enumerate(colors):
            for col, color in enumerate(values):
                c.rect(x + col * 190 / 7, 170 + row * 29, 190 / 7, 29, color,
                       stroke="#eaf2fa", sw=.45, id=f"{name}-row-{row + 1}-column-{col + 1}")


def inputs(c):
    with c.group("input-data", "Single-cell multi-omics inputs", layer=True):
        c.rect(20, 19, 482, 738, "url(#input-panel)", 20, "#b9bcc8", 1.8)
        with c.group("input-content", "Centered input contents"):
            c.parent.set("transform", "translate(0,24)")
            c.rect(31.5, 35, 459, 66, "url(#caption-fill)", 20)
            c.text(261, 79, "Single-cell multi-omics", 40, "bold")
            c.text(148.5, 157, "RNA", 33, "bold")
            c.text(373.5, 157, "ATAC", 33, "bold")
            matrix(c, "rna", 53.5, RNA_COLORS)
            matrix(c, "atac", 278.5, ATAC_COLORS)
            c.text(148.5, 321, "Genes × cells", 28)
            c.text(373.5, 321, "Regions × cells", 28)
            with c.group("paired-or-aligned", "Paired or aligned profiles"):
                c.rect(101, 348, 320, 49, "url(#caption-fill)", 21)
                c.text(261, 381, "Paired or aligned", 29)
            with c.group("candidate-resources", "Motifs and annotation candidate resources"):
                c.parent.set("transform", "translate(42,12)")
                c.text(219, 496, "Motifs + annotations", 35)
                with c.group("motif-logo", "Editable motif letters"):
                    c.parent.set("transform", "translate(28,14)")
                    for letter, x, y, size, color, length in [
                        ("T", 56, 596, 110, "#ff0711", 34),
                        ("G", 85, 572, 61, "#ffad00", 29),
                        ("A", 89, 596, 35, "#00b840", 26),
                        ("T", 116, 596, 110, "#ff0711", 33),
                        ("T", 150, 576, 65, "#ff0711", 27),
                        ("G", 151, 596, 34, "#ffad00", 25),
                        ("T", 179, 576, 52, "#ff0711", 23),
                        ("G", 181, 596, 34, "#ffad00", 24),
                    ]:
                        c.text(x, y, letter, size, "bold", "start", color, length=length)
                with c.group("annotation-database", "Annotation database cylinder"):
                    c.parent.set("transform", "translate(-4,2)")
                    c.path("M277,531 V604 C277,627 360,627 360,604 V531", INK, 4.3, "#ffffff")
                    c.path("M277,555 C277,578 360,578 360,555", INK, 4)
                    c.path("M277,581 C277,604 360,604 360,581", INK, 4)
                    c.element("ellipse", cx=318.5, cy=531, rx=41.5, ry=17,
                              fill="#ffffff", stroke=INK, stroke_width=4.3)


def method(c):
    with c.group("inference-framework", "Independent fits and region matching", layer=True):
        c.rect(550, 16, 704, 743, "url(#method-panel)", 19)
        with c.group("method-content", "Inference contents"):
            c.parent.set("transform", "translate(0,24)")
            c.text(909, 67, "Independent sparse fits", 39, "bold")
            math_text = c.text(902, 118, "Greedy–ℓ₀ · ", 31)
            ET.SubElement(math_text, tag("tspan"), {"font-style": "italic"}).text = "BIC"
            with c.group("network-flow", "Regulatory fits and application connections"):
                c.parent.set("transform", "translate(0,28)")
                with c.group("independent-tf-gene-fit", "Independent TF–gene fit"):
                    c.rect(567, 160, 201, 121, "#ffffff", 19)
                    c.text(664, 202, "TF–gene", 29)
                    c.tf(612, 240, 42, 40)
                    c.arrow(644, 241, 697, 241, head=18)
                    c.gene(724, 241, 23)
                with c.group("independent-tf-region-fit", "Independent TF–region fit"):
                    c.rect(787, 160, 220, 121, "#ffffff", 19)
                    c.text(899, 202, "TF–region", 29)
                    c.tf(834, 240, 41, 40)
                    c.arrow(866, 241, 921, 241, head=18)
                    c.diamond(953, 241, 25, 24)
                with c.group("independent-region-gene-fit", "Independent Region–gene fit"):
                    c.rect(1025, 160, 207, 121, "#ffffff", 19)
                    c.text(1128, 202, "Region–gene", 29)
                    c.diamond(1070, 241, 25, 24)
                    c.arrow(1104, 241, 1160, 241, head=18)
                    c.gene(1188, 241, 23)
                with c.group("exact-region-matching", "Only two chromatin layers enter exact matching"):
                    c.path("M875,282 V315 Q875,334 894,334 H921 Q941,334 941,355")
                    c.path("M1109,282 V315 Q1109,334 1090,334 H1063 Q1043,334 1043,355")
                    c.rect(807, 357, 309, 48, "url(#caption-fill)", 22)
                    c.text(962, 391, "Exact region matching", 28)
                    c.arrow(943, 407, 943, 455, head=20)
                with c.group("region-resolved-chain", "Region-resolved TF→Region→Gene chain"):
                    c.tf(743, 484, 56, 50)
                    c.text(742, 543, "TF", 27, "bold")
                    c.arrow(783, 484, 903, 484, head=22)
                    c.diamond(943, 484, 31, 29)
                    c.text(943, 540, "Region", 29, "bold")
                    c.arrow(986, 484, 1123, 484, head=22)
                    c.gene(1161, 484, 28.5)
                    c.text(1158, 545, "Gene", 28, "bold")
                    with c.group("genomic-region-strip", "Genomic region peaks"):
                        for index, d in enumerate([
                            "M861,586 L864,573 L868,571 L869,559 L883,559 L885,573 L889,586 Z",
                            "M928,586 L930,572 L933,571 L934,559 L946,559 L948,573 L951,586 Z",
                            "M972,586 L975,575 L994,575 L997,586 Z",
                            "M1000,586 L1003,573 L1007,571 L1008,559 L1020,559 L1022,574 L1025,586 Z",
                        ]):
                            c.path(d, None, fill="url(#region-node)", id=f"genomic-peak-{index + 1}")
                        c.path("M824,588 H1062", EDGE, 3.2)
                        c.text(945, 619, "Genomic region", 26)
                with c.group("direct-network-output", "Direct TF–gene output to erythroid"):
                    c.path("M657,282 V553")
                    c.rect(559, 555, 193, 75, "url(#caption-fill)", 18)
                    c.text(655, 590, "Direct TF–gene", 27)
                    c.text(655, 618, "networks", 27)
                    c.path("M657,631 V653 Q657,663 667,663 H1254")
                    c.arrow(1254, 663, 1282, 663, width=4.3, head=22)
                with c.group("region-chain-output", "Region chain output to fetal cortex"):
                    c.path("M1198,484 H1239 Q1249,484 1249,474 V171")
                with c.group("direct-network-cortex-output", "Direct TF–gene overlap comparison in fetal cortex"):
                    c.path("M768,220 H778 V150 Q778,140 788,140 H1239 Q1249,140 1249,150 V171")
                    c.text(1085, 132, "Direct-link overlap", 21)
                with c.group("cortex-network-output", "Both network outputs enter the fetal-cortex analysis"):
                    c.path("M1249,171 H1254")
                    c.arrow(1254, 171, 1282, 171, width=4.3, head=22)


NEURON_BRANCHES = {
    "blue": [
        ([(1386,96),(1384,108),(1380,115)],1.4),
        ([(1345,112),(1350,122)],1.4),
        ([(1380,116),(1381,129),(1382,130)],2.2),
        ([(1403,117),(1392,126),(1383,130)],2.2),
        ([(1351,123),(1351,125),(1355,132),(1359,136),(1359,140),(1363,148)],1.8),
        ([(1412,130),(1405,141),(1405,144)],2.2),
        ([(1382,131),(1382,134),(1379,138),(1378,143)],3.5),
        ([(1344,141),(1346,147)],1.4),
        ([(1422,142),(1411,145),(1406,145)],1.4),
        ([(1377,144),(1373,144),(1364,148)],1.4),
        ([(1379,144),(1384,146),(1391,152)],1.4),
        ([(1405,145),(1402,148),(1392,152)],2.2),
        ([(1347,148),(1360,152)],1.4),
        ([(1346,149),(1339,153),(1330,152),(1322,155),(1317,155)],1.4),
        ([(1363,149),(1361,152)],1.4),
        ([(1361,153),(1360,154),(1360,165),(1362,171)],2.2),
        ([(1391,153),(1391,158),(1392,159),(1391,167)],1.4),
        ([(1390,168),(1382,175)],1.4),
        ([(1392,168),(1397,169),(1400,172),(1408,173),(1417,176)],1.4),
        ([(1361,172),(1351,180),(1328,203)],3.8),
        ([(1363,172),(1375,176),(1380,175)],1.4),
        ([(1381,176),(1383,189),(1389,196),(1392,203),(1394,219)],1.4),
        ([(1418,177),(1427,187)],2.2),
        ([(1399,183),(1403,188)],1.4),
        ([(1307,192),(1323,199),(1326,202)],1.8),
        ([(1326,203),(1323,206),(1316,209)],2.2),
        ([(1328,203),(1330,205),(1330,207),(1329,208),(1329,216),(1327,225)],2.2),
        ([(1314,209),(1306,210)],1.4),
        ([(1315,210),(1314,215),(1310,219)],2.2),
    ],
    "green": [
        ([(1724,91),(1722,93),(1722,102),(1721,106),(1719,108)],1.4),
        ([(1715,103),(1718,108)],1.4),
        ([(1682,109),(1687,118)],1.4),
        ([(1719,109),(1719,125),(1718,126),(1719,130)],1.4),
        ([(1692,114),(1688,118)],1.4),
        ([(1742,115),(1725,128),(1719,130)],1.4),
        ([(1688,119),(1688,122),(1696,135),(1698,142),(1701,146)],1.4),
        ([(1749,129),(1743,139),(1743,143)],1.4),
        ([(1719,130),(1717,132),(1717,141)],3.5),
        ([(1679,139),(1682,146)],1.4),
        ([(1762,141),(1756,141),(1744,144)],1.4),
        ([(1716,142),(1708,143),(1706,145),(1702,146)],3.5),
        ([(1718,142),(1723,144),(1729,151)],3.8),
        ([(1742,144),(1739,147),(1730,151)],2.2),
        ([(1682,146),(1686,146),(1699,151)],1.4),
        ([(1682,146),(1678,150),(1664,150),(1655,152)],1.4),
        ([(1701,147),(1699,151)],3.5),
        ([(1699,151),(1698,153),(1698,162),(1701,169)],1.4),
        ([(1730,152),(1729,165)],2.8),
        ([(1729,166),(1737,168),(1740,171),(1757,174)],1.4),
        ([(1728,167),(1725,171),(1720,174)],1.4),
        ([(1702,170),(1709,174),(1713,175),(1718,174)],1.4),
        ([(1701,171),(1681,187),(1666,202)],3.8),
        ([(1719,175),(1721,187),(1725,191),(1730,202),(1731,218),(1732,219)],1.4),
        ([(1758,175),(1765,184),(1768,186)],1.4),
        ([(1738,183),(1742,187)],1.4),
        ([(1642,192),(1652,194),(1662,199),(1665,202)],1.4),
        ([(1664,203),(1654,209)],2.2),
        ([(1666,203),(1667,204),(1667,213),(1666,214),(1665,224)],1.4),
        ([(1652,209),(1641,210)],1.4),
        ([(1653,210),(1652,215),(1648,219)],2.2),
    ],
}


def neuron(c, name, x, y, green=False):
    color = "#00b74f" if green else "#009bfa"
    branches = NEURON_BRANCHES["green" if green else "blue"]
    origin = (1714, 159) if green else (1376, 160)
    with c.group(name, "Inhibitory neuron" if green else "Excitatory neuron"):
        for index, (points, width) in enumerate(branches):
            d = "M" + " L".join(f"{x + px - origin[0]:g},{y + py - origin[1]:g}" for px, py in points)
            c.path(d, color, width, id=f"{name}-branch-{index + 1}")
        c.circle(x, y, 17.2, "url(#neuron-green)" if green else "url(#neuron-blue)", color, 1.1)
        c.circle(x, y + .5, 8.6, "#2d9b75" if green else "#088ad6", sw=0)
        c.circle(x - 3, y - 3, 3.5, "#67d197" if green else "#48c6f5", opacity=.4)


def cortex_network(c, name, shift=0):
    with c.group(name, "TF–region–gene cortical example"):
        tf_positions = [1820, 1868, 1917] if shift else [1473.5, 1524.5, 1578]
        region_positions = [1829, 1903] if shift else [1487.5, 1567]
        gene_position = 1866 if shift else 1527.5
        for x in tf_positions:
            c.tf(x, 111, 28, 29)
        for x in region_positions:
            c.diamond(x, 185, 21, 20)
        c.gene(gene_position, 271, 22)
        left_tf, center_tf, right_tf = tf_positions
        left_region, right_region = region_positions
        for index, (x1, y1, x2, y2) in enumerate([
            (left_tf + 2, 129, left_region - 2, 164),
            (center_tf - 7, 130, left_region + 4, 162),
            (center_tf + 10, 130, right_region - 4, 163),
            (right_tf - 2, 129, right_region + 2, 164),
            (left_region + 8, 209, gene_position - 8, 249),
            (right_region - 8, 209, gene_position + 9, 249),
        ]):
            with c.group(f"{name}-edge-{index + 1}"):
                c.arrow(x1, y1, x2, y2, width=2.4, head=13)


def cortex(c):
    with c.group("fetal-cortex", "Fetal cortex examples", layer=True):
        c.rect(1283, 16, 681, 335, "url(#application-panel)", 18)
        c.text(1305, 56, "Fetal cortex", 28, "bold", "start")
        c.path("M1631,92 V331", "#b5b7c2", 1.2)
        neuron(c, "excitatory-neuron", 1376, 160)
        neuron(c, "inhibitory-neuron", 1714, 159, green=True)
        cortex_network(c, "znf804a-network")
        cortex_network(c, "kcnc2-network", 340)
        c.text(1371, 270, "Excitatory", 27)
        c.text(1376, 298, "neuron", 27)
        c.text(1717, 270, "Inhibitory", 27)
        c.text(1720, 299, "neuron", 27)
        c.text(1531, 326, "ZNF804A", 26, italic=True)
        c.text(1869, 326, "KCNC2", 26, italic=True)


def erythroid_cell(c, name, x, y, radius, nucleus):
    with c.group(name, "Nucleated erythroid cell"):
        c.circle(x, y, radius, "url(#cell-body)", "#ff1825", 2.1)
        c.circle(x, y, radius - 2.4, "none", "#fff6fc", 1.2)
        c.circle(x, y, nucleus, "url(#cell-nucleus)")


def erythroid(c):
    with c.group("erythroid-maturation", "Erythroid direct-network target turnover", layer=True):
        c.rect(1283, 365, 681, 393, "url(#application-panel)", 18)
        c.text(1303, 404, "Erythroid maturation", 31, "bold", "start")
        c.path("M1597,447 V689", "#b5b7c2", 1.2)
        c.path("M1454,470 V604", "#b5b7c2", 1.1)
        c.path("M1854,470 V604", "#b5b7c2", 1.1)
        with c.group("hmga1-example", "HMGA1 Erythroblast-to-Normoblast lost targets"):
            erythroid_cell(c, "hmga1-erythroblast", 1335, 536, 28, 13)
            erythroid_cell(c, "hmga1-normoblast", 1419, 538, 25.5, 11.5)
            c.arrow(1368, 536, 1389, 536, width=2.6, head=11)
            c.text(1379, 591, "Ery → Norm", 22)
            c.text(1527, 470, "HMGA1", 26, italic=True)
            c.tf(1522, 493, 28, 27)
            c.arrow(1515, 511, 1498, 557, width=2.4, head=13, dash="9 7")
            c.arrow(1532, 511, 1550, 557, width=2.4, head=13, dash="9 7")
            c.gene(1498, 575, 17.5)
            c.gene(1550, 575, 17.5)
            c.text(1443, 648, "Lost targets", 24, "bold")
            c.text(1440, 678, "Translation / mitochondria", 24)
        with c.group("bnc2-example", "BNC2 targets gained and lost across successive transitions"):
            erythroid_cell(c, "bnc2-proerythroblast", 1642, 537, 26, 12)
            erythroid_cell(c, "bnc2-erythroblast", 1726, 537, 25, 12)
            erythroid_cell(c, "bnc2-normoblast", 1811, 537, 26, 11.5)
            c.arrow(1673, 537, 1697, 537, width=2.5, head=11)
            c.arrow(1757, 537, 1783, 537, width=2.5, head=11)
            c.text(1684, 496, "Gained", 21, color="#008b16")
            c.path("M1664,502 L1657,513", "#4aa850", 1.4)
            c.path("M1706,502 L1714,513", "#4aa850", 1.4)
            c.text(1795, 496, "Lost", 21)
            for x, label in [(1642, "ProE"), (1726, "Ery"), (1811, "Norm")]:
                c.text(x, 590, label, 22)
            c.text(1907, 470, "BNC2", 26, italic=True)
            c.tf(1907, 494, 28, 27)
            c.arrow(1900, 515, 1883, 557, GREEN, 2.4, 13)
            c.arrow(1918, 515, 1936, 557, width=2.4, head=13, dash="9 7")
            c.gene(1883, 575, 17.5)
            c.gene(1936, 575, 17.5)
            c.text(1779, 648, "Gained / lost targets", 24, "bold")
            c.text(1785, 678, "Immune responses", 24)
        with c.group("edge-style-key", "Gained and lost relation styles"):
            c.arrow(1473, 725, 1535, 725, GREEN, 2.5, 16)
            c.text(1552, 732, "Gained", 21, anchor="start")
            c.arrow(1660, 725, 1722, 725, EDGE, 2.5, 16, "10 7")
            c.text(1739, 732, "Lost", 21, anchor="start")


def draw():
    c = Canvas()
    definitions(c)
    c.rect(0, 0, WIDTH, HEIGHT, "#ffffff", id="canvas-background")
    inputs(c)
    method(c)
    cortex(c)
    erythroid(c)
    with c.group("common-input-arrow", "All inputs enter the inference framework", layer=True):
        c.arrow(502, 351, 549, 351, width=4.3, head=22)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    ET.indent(c.root, space="  ")
    ET.ElementTree(c.root).write(OUTPUT, encoding="utf-8", xml_declaration=True)
    for format_name in ("png", "pdf"):
        subprocess.run([
            "rsvg-convert", "--format", format_name,
            "--output", str(OUTPUT.with_suffix(f".{format_name}")), str(OUTPUT),
        ], check=True)
    return OUTPUT


if __name__ == "__main__":
    print(draw())
