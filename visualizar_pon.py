import json, math, sys
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

COLORES_TASA = {100: "#9ecae1", 500: "#3182bd", 1000: "#08306b"}
COLORES_PUERTO = {"GPON": "#2ca02c", "XG-PON": "#ff7f0e", "XGS-PON": "#d62728"}
JERARQUIA_PUERTO = {"GPON": 0, "XG-PON": 1, "XGS-PON": 2}
ABREVIATURA = {"GPON": "G", "XG-PON": "XG", "XGS-PON": "XGS"}

def segmentos_de(claves, xy_nodos):
    segmentos = []
    for clave in claves:
        extremos = clave.split("|")
        if extremos[0] in xy_nodos and extremos[1] in xy_nodos:
            segmentos.append([xy_nodos[extremos[0]], xy_nodos[extremos[1]]])
    return segmentos

def graficar(mundo, resultado, sitio, semilla, ruta_png=None, mostrar=True):
    xy_nodos = mundo["xy_nodos"]
    xy_puntos = mundo["xy_puntos"]
    solucion = resultado["solucion"]
    figura, eje = plt.subplots(figsize=(13, 10))

    eje.add_collection(LineCollection(segmentos_de(mundo["largo"].keys(), xy_nodos),
                                      colors="#d9d9d9", linewidths=0.4, zorder=1))
    eje.add_collection(LineCollection(segmentos_de(solucion["zanja"], xy_nodos),
                                      colors="#444444", linewidths=2.2, zorder=3))

    enlaces = []
    for i in solucion["asignacion"]:
        destino = solucion["asignacion"][i]
        enlaces.append([xy_puntos[i], xy_nodos[destino[0]]])
    eje.add_collection(LineCollection(enlaces, colors="#e6550d", linewidths=0.8,
                                      linestyles="dashed", zorder=4))

    for tasa in COLORES_TASA:
        xs, ys, tamanos = [], [], []
        for i in mundo["n"]:
            if i in solucion["sin_cobertura"]:
                continue
            if solucion["tasa"][i] != tasa:
                continue
            xs.append(xy_puntos[i][0])
            ys.append(xy_puntos[i][1])
            tamanos.append(20 + 8 * math.sqrt(mundo["n"][i]))
        eje.scatter(xs, ys, s=tamanos, c=COLORES_TASA[tasa], edgecolors="black",
                    linewidths=0.5, zorder=5)

    xs_sin, ys_sin = [], []
    for i in solucion["sin_cobertura"]:
        xs_sin.append(xy_puntos[i][0])
        ys_sin.append(xy_puntos[i][1])
    if len(xs_sin) > 0:
        eje.scatter(xs_sin, ys_sin, marker="X", s=90, c="red", edgecolors="black",
                    linewidths=0.6, zorder=6)

    xs_cand, ys_cand = [], []
    for candidato in mundo["sitios"]:
        if candidato == sitio:
            continue
        xs_cand.append(xy_nodos[candidato][0])
        ys_cand.append(xy_nodos[candidato][1])
    eje.scatter(xs_cand, ys_cand, marker="^", s=45, facecolors="none",
                edgecolors="#555555", linewidths=0.9, zorder=6)

    por_nodo = {}
    for fila in solucion["splitters"]:
        por_nodo.setdefault(fila[0], []).append(fila)
    for j in por_nodo:
        mayor = "GPON"
        lineas = []
        for fila in por_nodo[j]:
            if JERARQUIA_PUERTO[fila[1]] > JERARQUIA_PUERTO[mayor]:
                mayor = fila[1]
            lineas.append(f"{ABREVIATURA[fila[1]]} 1:{fila[2]} x{fila[3]}")
        eje.scatter([xy_nodos[j][0]], [xy_nodos[j][1]], marker="s", s=110,
                    c=COLORES_PUERTO[mayor], edgecolors="black", linewidths=1.0, zorder=8)
        eje.annotate("\n".join(lineas), xy_nodos[j], xytext=(6, 6), textcoords="offset points",
                     fontsize=6.5, zorder=9,
                     bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))

    eje.scatter([xy_nodos[sitio][0]], [xy_nodos[sitio][1]], marker="*", s=600,
                c="gold", edgecolors="black", linewidths=1.2, zorder=10)

    todos_x, todos_y = [], []
    for i in xy_puntos:
        todos_x.append(xy_puntos[i][0])
        todos_y.append(xy_puntos[i][1])
    margen = 1.5
    eje.set_xlim(min(todos_x) - margen, max(todos_x) + margen)
    eje.set_ylim(min(todos_y) - margen, max(todos_y) + margen)
    eje.set_aspect("equal")
    eje.set_xlabel("Este (km)")
    eje.set_ylabel("Norte (km)")

    total_clientes = sum(mundo["n"].values())
    porcentaje = 100 * resultado["cob"] / total_clientes
    eje.set_title(f"Semilla {semilla} | OLT en {sitio} | cobertura {porcentaje:.1f}% | "
                  f"costo {resultado['costo']:.1f} | puertos {resultado['puertos']} | {resultado['estado']}")

    leyenda = [
        Line2D([0], [0], marker="*", color="w", markerfacecolor="gold", markeredgecolor="black",
               markersize=18, label="OLT elegida"),
        Line2D([0], [0], marker="^", color="w", markerfacecolor="none", markeredgecolor="#555555",
               markersize=9, label="sitio candidato"),
        Line2D([0], [0], color="#444444", linewidth=2.2, label="zanja"),
        Line2D([0], [0], color="#e6550d", linewidth=1, linestyle="dashed", label="fibra de distribución"),
        Line2D([0], [0], marker="X", color="w", markerfacecolor="red", markeredgecolor="black",
               markersize=10, label="sin cobertura"),
    ]
    for tasa in COLORES_TASA:
        leyenda.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=COLORES_TASA[tasa],
                              markeredgecolor="black", markersize=9, label=f"demanda {tasa} Mbps"))
    for tipo in COLORES_PUERTO:
        leyenda.append(Line2D([0], [0], marker="s", color="w", markerfacecolor=COLORES_PUERTO[tipo],
                              markeredgecolor="black", markersize=9, label=f"splitter {tipo}"))
    eje.legend(handles=leyenda, loc="upper left", fontsize=8, framealpha=0.9)

    figura.tight_layout()
    if ruta_png is not None:
        figura.savefig(ruta_png, dpi=200)
    if mostrar:
        plt.show()
    plt.close(figura)

if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8") as archivo:
        mundo = json.load(archivo)
    with open(sys.argv[2], encoding="utf-8") as archivo:
        guardado = json.load(archivo)
    graficar(mundo, guardado["resultado"], guardado["sitio"], guardado["semilla"],
             f"red_pon_semilla{guardado['semilla']}.png", mostrar=True)