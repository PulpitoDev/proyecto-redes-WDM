# === Librerías ===
import json, math, random, time, sys
import networkx as nx, pulp as pl

# === Parámetros en CLP (fuentes citadas; ver supuestos en el informe) ===
N_TOTAL = 3000                                   # clientes totales (caso base)
NIVELES = [100, 500, 1000]                       # tasa mínima simétrica por zona (Mbps)
ESCENARIO_PRECIO = "medio"                       # "min", "medio" o "max": rango de precio de las tarjetas de puertos
FILA = 16                                        # puertos por tarjeta (tarjetas de 16 puertos)
PRECIO_TARJETA = {                               # CLP por tarjeta de 16 puertos: (min, max)  [router-switch.com, Huawei SmartAX MA5800]
    "GPON":    (2_500_000, 2_500_000),           # H901GPHF
    "XG-PON":  (2_500_000, 4_500_000),           # H901XGSF
    "XGS-PON": (3_000_000, 5_500_000)}           # H901XSHF
TASAS = {"GPON": (2488, 1244), "XG-PON": (9953, 2488), "XGS-PON": (9953, 9953)}   # Mbps de línea (bajada, subida)
def precio_tarjeta(t):
    lo, hi = PRECIO_TARJETA[t]
    return {"min": lo, "max": hi}.get(ESCENARIO_PRECIO, (lo + hi) / 2)
PUERTOS = {t: (TASAS[t][0], TASAS[t][1]) for t in TASAS}
TARJETAS_MAX = 17                                # 17 tarjetas x 16 puertos = 272 puertos (chasis MA5800-X17)
K_R = {8: 22_610.0, 16: 35_105.0}                # CLP por splitter 1:r  [dato del grupo]
LMAX = {8: 20.0, 16: 20.0}                       # alcance por razón (km), presupuesto de potencia
DMAX = 20.0                                      # alcance máximo (km)
BETA, ALFA = 2_379_000.0, 0                     # CLP por km: BETA = cable de fibra [compratecno.cl], ALFA = zanja (estimada)
UNIDAD = 1e6                                     # el solver trabaja en millones de CLP
CONFIGS = [(t, r) for t in PUERTOS for r in K_R]
LIMITE_S = 60                                    # segundos máximos por MILP

# === Estructura "mundo" (la misma para datos sintéticos y reales) ===
# n[i], zona[i]: clientes y zona de cada punto de demanda; J: nodos candidatos de splitter; sitios: candidatos de OLT
# dsj[s][j]: km por calle OLT->nodo ; dji[j][i]: km por calle nodo->punto ; rutas[s][j]: tramos de la ruta ; largo[e]: km del tramo

def mundo_sintetico(n_puntos=24):
    rng = random.Random(0)
    G = nx.grid_2d_graph(10, 8)
    for a, b in G.edges: G[a][b]["l"] = round(rng.uniform(1.0, 2.0), 2)
    riesgo = {v for v in G if v[0] >= 8 and v[1] <= 2}
    puntos = rng.sample(sorted(G.nodes), n_puntos)
    w = [rng.uniform(0.5, 1.5) for _ in puntos]
    n = [max(1, round(N_TOTAL * x / sum(w))) for x in w]
    n[0] += N_TOTAL - sum(n)
    nm = lambda v: f"{v[0]}_{v[1]}"
    ek = lambda a, b: "|".join(sorted((nm(a), nm(b))))
    J = [v for v in G if v not in riesgo and (v[0] + v[1]) % 2 == 0]
    sitios = [v for v in G if v[0] % 2 == 0 and v[1] % 2 == 0 and v not in riesgo]
    dist = dict(nx.all_pairs_dijkstra_path_length(G, weight="l"))
    mundo = dict(n={nm(p): q for p, q in zip(puntos, n)}, zona={nm(p): f"{p[0]//3}_{p[1]//3}" for p in puntos},
                 J=[nm(j) for j in J], sitios=[nm(s) for s in sitios], dsj={}, dji={}, rutas={},
                 largo={ek(a, b): G[a][b]["l"] for a, b in G.edges})
    for s in sitios:
        paths = nx.single_source_dijkstra_path(G, s, weight="l")
        mundo["dsj"][nm(s)] = {nm(j): dist[s][j] for j in J}
        mundo["rutas"][nm(s)] = {nm(j): [ek(a, b) for a, b in zip(paths[j], paths[j][1:])] for j in J}
    for j in J:
        mundo["dji"][nm(j)] = {nm(p): dist[j][p] for p in puntos if dist[j][p] <= 3.5}
    return mundo

def cargar_mundo(ruta):
    with open(ruta, encoding="utf-8") as f: return json.load(f)

# === MILP de diseño de la PON para una OLT fija en s ===
def resolver(s, mundo, m, limite=LIMITE_S):
    n, J, dsj, dji, rutas, largo = (mundo[k] for k in ("n", "J", "dsj", "dji", "rutas", "largo"))
    pares = {}                                                          # (i,j) -> configuraciones admisibles
    for j in J:
        if j not in dsj[s]: continue
        for i, dij in dji[j].items():
            cs = [(t, r) for (t, r) in CONFIGS
                  if PUERTOS[t][0] / r >= m[i] and PUERTOS[t][1] / r >= m[i]
                  and dsj[s][j] + dij <= min(DMAX, LMAX[r])]
            if cs: pares[(i, j)] = cs
    cubiertos = {i for (i, j) in pares}
    cob = sum(n[i] for i in cubiertos)
    if not cubiertos: return dict(estado="sin cobertura", cob=0, costo=None)
    Ju = sorted({j for (_, j) in pares})
    tern = [(i, j, c) for (i, j), cs in pares.items() for c in cs]
    V = lambda nom, idx, **kw: {k: pl.LpVariable(f"{nom}{q}", **kw) for q, k in enumerate(idx)}
    a = V("a", tern, cat="Binary")
    x = V("x", [(j, c) for j in Ju for c in CONFIGS], lowBound=0, cat="Integer")
    w = V("w", Ju, cat="Binary")
    E = sorted({e for j in Ju for e in rutas[s][j]})
    g = V("g", E, cat="Binary")
    mod = pl.LpProblem("pon", pl.LpMinimize)
    dji_ij = {(i, j): dji[j][i] for (i, j) in pares}
    tar = V("tar", list(PUERTOS), lowBound=0, cat="Integer")                                  # tarjetas de 16 puertos por tecnología
    mod += (pl.lpSum(K_R[r] / UNIDAD * x[(j, (t, r))] for j in Ju for (t, r) in CONFIGS)       # splitters
            + pl.lpSum(precio_tarjeta(t) / UNIDAD * tar[t] for t in PUERTOS)                    # tarjetas de la OLT
            + pl.lpSum((ALFA + BETA) / UNIDAD * largo[e] * g[e] for e in E)                     # alimentación: zanja + un cable por tramo
            + pl.lpSum((ALFA + BETA) / UNIDAD * dji_ij[(i, j)] * a[(i, j, c)] for (i, j, c) in tern))   # distribución: zanja + cable por punto
    por_i, por_jc = {}, {}
    for k in tern:
        por_i.setdefault(k[0], []).append(k); por_jc.setdefault((k[1], k[2]), []).append(k)
    for i in cubiertos: mod += pl.lpSum(a[k] for k in por_i[i]) == 1                          # (1) asignación única
    for (j, c), ks in por_jc.items():
        mod += pl.lpSum(n[k[0]] * a[k] for k in ks) <= c[1] * x[(j, c)]                      # (2) capacidad
        for k in ks: mod += a[k] <= w[j]                                                      # (3) nodo activo
    for j in Ju:
        for e in rutas[s][j]: mod += g[e] >= w[j]                                             # (4) zanja compartida
    for t in PUERTOS:                                                                         # (5) puertos <= 16 por tarjeta
        mod += pl.lpSum(x[(j, c)] for j in Ju for c in CONFIGS if c[0] == t) <= FILA * tar[t]
    mod += pl.lpSum(tar.values()) <= TARJETAS_MAX                                             # (6) ranuras de la OLT
    t0 = time.time()
    mod.solve(pl.PULP_CBC_CMD(msg=False, timeLimit=limite, gapRel=0.02))
    est = pl.LpStatus[mod.status]
    if est != "Optimal": return dict(estado=est, cob=cob, costo=None, t=time.time() - t0)
    est = "Optimo" if mod.sol_status == 1 else "Factible"          # Factible = límite de tiempo, sin prueba de optimalidad
    sirve = {k[0]: k for k in tern if a[k].value() > 0.5}
    assert set(sirve) == cubiertos, "punto sin servir"
    carga = {}
    for i, (_, j, c) in sirve.items():
        assert c in pares[(i, j)], "configuración no admisible (QoS o alcance)"
        carga[(j, c)] = carga.get((j, c), 0) + n[i]
    for (j, c), q in carga.items(): assert q <= c[1] * round(x[(j, c)].value()), "capacidad violada"
    puertos = sum(round(v.value()) for v in x.values())
    assert puertos <= FILA * TARJETAS_MAX
    sp = {}
    for (j, c), v in x.items():
        if round(v.value()) > 0: sp.setdefault(j, []).append(dict(tipo=c[0], razon=c[1], cantidad=round(v.value())))
    diseno = dict(splitters=sp, asignacion={i: dict(nodo=j, tipo=c[0], razon=c[1]) for i, (_, j, c) in sirve.items()},
                  tramos=[e for e in E if g[e].value() > 0.5])
    cuenta = lambda f: sum(round(x[k].value()) for k in x if f(k[1]))
    return dict(estado=est, cob=cob, costo=pl.value(mod.objective) * UNIDAD, puertos=puertos, tarjetas=round(sum(v.value() for v in tar.values())), t=time.time() - t0,
                gpon=cuenta(lambda c: c[0] == "GPON"), xg=cuenta(lambda c: c[0] == "XG-PON"),
                xgs=cuenta(lambda c: c[0] == "XGS-PON"), r8=cuenta(lambda c: c[1] == 8), diseno=diseno)

# === Barrido de sitios para una semilla ===
TOL_COB = 0.005                                   # tolerancia de cobertura (0,5 puntos)

def barrido(semilla, mundo, verbose=True):
    rng = random.Random(semilla)
    zonas = sorted(set(mundo["zona"].values()))
    mz = {z: rng.choice(NIVELES) for z in zonas}                       # tasa mínima por zona (sorteo con semilla)
    m = {i: mz[z] for i, z in mundo["zona"].items()}
    res = {s: resolver(s, mundo, m) for s in mundo["sitios"]}
    ok = {s: r for s, r in res.items() if r["costo"] is not None}
    mejor = None
    if ok:                                                             # cobertura con tolerancia, luego menor costo
        mejor_cob = max(r["cob"] for r in ok.values())
        cand = {s: r for s, r in ok.items() if r["cob"] >= mejor_cob - TOL_COB * N_TOTAL}
        mejor = min(cand, key=lambda s: cand[s]["costo"])
    if verbose:
        print(f"Semilla {semilla}: zonas con 100/500/1000 Mbps =", [list(mz.values()).count(v) for v in NIVELES])
        for s in mundo["sitios"]:
            r = res[s]
            print(f"  sitio {s}: cobertura {100 * r['cob'] / N_TOTAL:5.1f}% | estado {r['estado']:10s} |",
                  "costo ---" if r["costo"] is None else f"costo {r['costo'] / 1e6:8.1f} M$ | puertos {r['puertos']:3d} | tarjetas {r['tarjetas']:2d} | t {r['t']:.1f}s")
        print("  -> mejor sitio:", mejor)
    return mejor, (ok.get(mejor) if mejor else None)

if __name__ == "__main__":
    # uso: python modelo_pon.py sintetico 1 2 3   |   python modelo_pon.py datos_pon.json 1 2 3
    fuente, semillas = sys.argv[1], [int(a) for a in sys.argv[2:]] or [1]
    mundo = mundo_sintetico() if fuente == "sintetico" else cargar_mundo(fuente)
    N_TOTAL = sum(mundo["n"].values())
    for sm in semillas: barrido(sm, mundo)
