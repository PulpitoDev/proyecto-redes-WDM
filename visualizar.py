# === Librerías ===
import sys, json, random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import modelo_pon as M

# === Uso: py visualizar.py [semilla] [sitio] ===  (sin sitio: usa el mejor del barrido)
semilla = int(sys.argv[1]) if len(sys.argv) > 1 else 1
sitio = sys.argv[2] if len(sys.argv) > 2 else None
mundo = M.cargar_mundo("datos_pon.json")
M.N_TOTAL = sum(mundo["n"].values())
assert "coords" in mundo, "datos_pon.json no tiene coordenadas: vuelve a correr datos_reales.py"

# === Tasas por zona (mismo sorteo que el barrido) ===
rng = random.Random(semilla)
mz = {z: rng.choice(M.NIVELES) for z in sorted(set(mundo["zona"].values()))}
m = {i: mz[z] for i, z in mundo["zona"].items()}
if sitio is None: sitio, _ = M.barrido(semilla, mundo, verbose=False)
r = M.resolver(sitio, mundo, m)
d, C, XY, PXY = r["diseno"], mundo["n"], mundo["coords"], mundo["puntos_xy"]
req = sorted({e for j in d["splitters"] for e in mundo["rutas"][sitio][j]})   # tramos realmente necesarios (excluye restos del solver)

# === Resumen del diseño ===
sp = d["splitters"]
total = sum(c["cantidad"] for L in sp.values() for c in L)
print(f"Semilla {semilla} | OLT en {sitio} | cobertura {100 * r['cob'] / M.N_TOTAL:.1f}% | costo {r['costo'] / 1e6:.1f} M$ CLP | {r['estado']}")
print(f"nodos con splitters: {len(sp)} | splitters (= puertos de OLT): {total}")
print("splitters por tipo y razón:")
for t in M.PUERTOS:
    for rz in sorted(M.K_R):
        q = sum(c["cantidad"] for L in sp.values() for c in L if c["tipo"] == t and c["razon"] == rz)
        if q: print(f"  {t:8s} 1:{rz:<2d} -> {q}")
print(f"zanja de alimentación: {sum(mundo['largo'][e] for e in req):.1f} km en {len(req)} tramos")
print(f"fibra de alimentación: {sum(M.BETA and mundo['dsj'][sitio][j] * sum(c['cantidad'] for c in L) for j, L in sp.items()):.0f} km-fibra | "
      f"fibra de distribución: {sum(C[i] * mundo['dji'][a['nodo']][i] for i, a in d['asignacion'].items()):.0f} km-fibra")
serv = {}
for i, a in d["asignacion"].items(): serv[a["nodo"]] = serv.get(a["nodo"], 0) + C[i]
print("nodos con más splitters (nodo | lon, lat | configuración | splitters | clientes):")
for j in sorted(sp, key=lambda j: -sum(c["cantidad"] for c in sp[j]))[:10]:
    conf = ", ".join(f"{c['cantidad']}x {c['tipo']} 1:{c['razon']}" for c in sp[j])
    print(f"  {j} | {XY[j][0]:.4f}, {XY[j][1]:.4f} | {conf} | {sum(c['cantidad'] for c in sp[j])} | {serv.get(j, 0)}")

# === Mapa ===
lat0 = np.mean([v[1] for v in PXY.values()])
fig, (a1, a2) = plt.subplots(1, 2, figsize=(19, 9.5))
segs = [(XY[u], XY[v]) for u, v in (e.split("|") for e in mundo["largo"])]
col = {"GPON": "#F2A541", "XG-PON": "#3E7CB1", "XGS-PON": "#00808A"}
for a in (a1, a2):
    a.add_collection(LineCollection(segs, colors="#C9D2DC", linewidths=0.3, zorder=1))
    a.set_aspect(1 / np.cos(np.radians(lat0))); a.set_xticks([]); a.set_yticks([])
xs, ys = [v[0] for v in PXY.values()], [v[1] for v in PXY.values()]
for a in (a1, a2):
    a.set_xlim(min(xs) - 0.01, max(xs) + 0.01); a.set_ylim(min(ys) - 0.01, max(ys) + 0.01)

# panel 1: zanja de alimentación, OLT y splitters
a1.add_collection(LineCollection([(XY[u], XY[v]) for u, v in (e.split("|") for e in req)], colors="#00808A", linewidths=2.2, zorder=3))
for j, L in sp.items():
    dom = max(L, key=lambda c: c["cantidad"])
    q = sum(c["cantidad"] for c in L)
    a1.scatter(*XY[j], s=25 + 7 * q, c=col[dom["tipo"]], edgecolors="#0B2545", linewidths=0.8, zorder=5)
    if q >= 8: a1.annotate(str(q), XY[j], ha="center", va="center", fontsize=7, color="white", zorder=6)
a1.scatter(*XY[sitio], marker="*", s=700, c="#C8553D", edgecolors="#0B2545", zorder=7)
a1.set_title(f"Alimentación: OLT, zanja y splitters (semilla {semilla})", fontsize=14, color="#0B2545")
a1.legend(handles=[Line2D([], [], marker="*", ls="", c="#C8553D", ms=16, label="OLT"),
                   Line2D([], [], c="#00808A", lw=2.5, label="zanja de alimentación")] +
                  [Line2D([], [], marker="o", ls="", c=c, mec="#0B2545", ms=9, label=f"splitters {t} (tamaño = cantidad)") for t, c in col.items()],
          loc="lower left", fontsize=9, framealpha=0.95)

# panel 2: clientes, tasa mínima y asignación
niv = {100: "#7FB77E", 500: "#F2A541", 1000: "#C8553D"}
for i, a in d["asignacion"].items():
    a2.plot([PXY[i][0], XY[a["nodo"]][0]], [PXY[i][1], XY[a["nodo"]][1]], c="#8DA9C4", lw=0.5, zorder=2)
for i in C:
    a2.scatter(*PXY[i], s=12 + C[i] * 1.2, c=niv[m[i]], edgecolors="#0B2545", linewidths=0.4, alpha=0.9, zorder=4)
for j in sp: a2.scatter(*XY[j], marker="^", s=45, c="#0B2545", zorder=5)
a2.scatter(*XY[sitio], marker="*", s=700, c="#C8553D", edgecolors="#0B2545", zorder=7)
a2.set_title("Distribución: clientes, tasa mínima y splitter asignado", fontsize=14, color="#0B2545")
a2.legend(handles=[Line2D([], [], marker="o", ls="", c=c, mec="#0B2545", ms=9, label=f"zona de {k} Mbps (tamaño = clientes)") for k, c in niv.items()] +
                  [Line2D([], [], marker="^", ls="", c="#0B2545", ms=8, label="nodo con splitters")], loc="lower left", fontsize=9, framealpha=0.95)
plt.tight_layout(); plt.savefig("diseno_pon.png", dpi=170, bbox_inches="tight")
json.dump(dict(semilla=semilla, sitio=sitio, costo=r["costo"], cobertura=r["cob"], diseno=d), open("diseno_pon.json", "w"))
print("Guardado: diseno_pon.png y diseno_pon.json")
plt.show()
