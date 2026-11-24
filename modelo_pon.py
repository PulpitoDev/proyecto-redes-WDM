import json, math, random, time, sys
import networkx as nx, pulp as pl
from visualizar_pon import graficar

n_total = 3000
NIVELES = [100, 500, 1000]
PUERTOS = {"GPON": (2488, 1244, 10.0),
           "XG-PON": (9953, 2488, 18.0),
           "XGS-PON": (9953, 9953, 25.0)}
K_R = {8: 0.8, 16: 1.2}
LMAX = {8: 20.0, 16: 14.0}
DMAX, P_MAX = 20.0, 400
BETA, ALFA = 0.5, 8.0
CONFIGS = [(t, r) for t in PUERTOS for r in K_R]
LIMITE_S = 10

def mundo_sintetico(n_puntos=24):
    rng = random.Random(0)
    G = nx.grid_2d_graph(10, 8)
    for a, b in G.edges: G[a][b]["l"] = round(rng.uniform(1.0, 2.0), 2)
    riesgo = {v for v in G if v[0] >= 8 and v[1] <= 2}
    puntos = rng.sample(sorted(G.nodes), n_puntos)
    w = [rng.uniform(0.5, 1.5) for _ in puntos]
    n = [max(1, round(n_total * x / sum(w))) for x in w]
    n[0] += n_total - sum(n)
    nm = lambda v: f"{v[0]}_{v[1]}"
    ek = lambda a, b: "|".join(sorted((nm(a), nm(b))))
    J = [v for v in G if v not in riesgo and (v[0] + v[1]) % 2 == 0]
    sitios = [v for v in G if v[0] % 2 == 0 and v[1] % 2 == 0 and v not in riesgo]
    dist = dict(nx.all_pairs_dijkstra_path_length(G, weight="l"))
    xy_nodos = {}
    for v in G.nodes:
        xy_nodos[nm(v)] = [float(v[0]), float(v[1])]
    xy_puntos = {}
    for p in puntos:
        xy_puntos[nm(p)] = [float(p[0]), float(p[1])]
    mundo = dict(n={nm(p): q for p, q in zip(puntos, n)}, zona={nm(p): f"{p[0]//3}_{p[1]//3}" for p in puntos},
                 J=[nm(j) for j in J], sitios=[nm(s) for s in sitios], dsj={}, dji={}, rutas={},
                 largo={ek(a, b): G[a][b]["l"] for a, b in G.edges},
                 xy_nodos=xy_nodos, xy_puntos=xy_puntos)
    for s in sitios:
        paths = nx.single_source_dijkstra_path(G, s, weight="l")
        mundo["dsj"][nm(s)] = {nm(j): dist[s][j] for j in J}
        mundo["rutas"][nm(s)] = {nm(j): [ek(a, b) for a, b in zip(paths[j], paths[j][1:])] for j in J}
    for j in J:
        mundo["dji"][nm(j)] = {nm(p): dist[j][p] for p in puntos if dist[j][p] <= 3.5}
    return mundo

def cargar_mundo(ruta):
    with open(ruta, encoding="utf-8") as f: return json.load(f)

def resolver(s, mundo, m, limite=LIMITE_S):
    n, J, dsj, dji, rutas, largo = (mundo[k] for k in ("n", "J", "dsj", "dji", "rutas", "largo"))
    pares = {}
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
    mod += (pl.lpSum((PUERTOS[t][2] + K_R[r] + BETA * dsj[s][j]) * x[(j, (t, r))] for j in Ju for (t, r) in CONFIGS)
            + pl.lpSum(ALFA * largo[e] * g[e] for e in E)
            + pl.lpSum(BETA * n[i] * dji_ij[(i, j)] * a[(i, j, c)] for (i, j, c) in tern))
    por_i, por_jc = {}, {}
    for k in tern:
        por_i.setdefault(k[0], []).append(k); por_jc.setdefault((k[1], k[2]), []).append(k)
    for i in cubiertos: mod += pl.lpSum(a[k] for k in por_i[i]) == 1
    for (j, c), ks in por_jc.items():
        mod += pl.lpSum(n[k[0]] * a[k] for k in ks) <= c[1] * x[(j, c)]
        for k in ks: mod += a[k] <= w[j]
    for j in Ju:
        for e in rutas[s][j]: mod += g[e] >= w[j]
    mod += pl.lpSum(x.values()) <= P_MAX
    t0 = time.time()
    mod.solve(pl.PULP_CBC_CMD(msg=False, timeLimit=limite, gapRel=0.02))
    est = pl.LpStatus[mod.status]
    if est != "Optimal": return dict(estado=est, cob=cob, costo=None, t=time.time() - t0)
    est = "Optimo" if mod.sol_status == 1 else "Factible"
    sirve = {k[0]: k for k in tern if a[k].value() > 0.5}
    assert set(sirve) == cubiertos, "punto sin servir"
    carga = {}
    for i, (_, j, c) in sirve.items():
        assert c in pares[(i, j)], "configuración no admisible (QoS o alcance)"
        carga[(j, c)] = carga.get((j, c), 0) + n[i]
    for (j, c), q in carga.items(): assert q <= c[1] * round(x[(j, c)].value()), "capacidad violada"
    puertos = sum(round(v.value()) for v in x.values())
    assert puertos <= P_MAX
    cuenta = lambda f: sum(round(x[k].value()) for k in x if f(k[1]))

    splitters = []
    for clave in x:
        cantidad = int(round(x[clave].value()))
        if cantidad > 0:
            splitters.append([clave[0], clave[1][0], clave[1][1], cantidad])
    asignacion = {}
    for i in sirve:
        terna = sirve[i]
        asignacion[i] = [terna[1], terna[2][0], terna[2][1]]
    zanja = []
    for e in E:
        if g[e].value() > 0.5:
            zanja.append(e)
    tasa = {}
    for i in n:
        tasa[i] = m[i]
    solucion = dict(sitio=s, splitters=splitters, asignacion=asignacion, zanja=zanja,
                    sin_cobertura=sorted(set(n) - cubiertos), tasa=tasa)

    return dict(estado=est, cob=cob, costo=pl.value(mod.objective), puertos=puertos, t=time.time() - t0,
                gpon=cuenta(lambda c: c[0] == "GPON"), xg=cuenta(lambda c: c[0] == "XG-PON"),
                xgs=cuenta(lambda c: c[0] == "XGS-PON"), r8=cuenta(lambda c: c[1] == 8),
                solucion=solucion)

TOL_COB = 0.005

def barrido(semilla, mundo, verbose=True):
    rng = random.Random(semilla)
    zonas = sorted(set(mundo["zona"].values()))
    mz = {z: rng.choice(NIVELES) for z in zonas}
    m = {i: mz[z] for i, z in mundo["zona"].items()}
    res = {s: resolver(s, mundo, m) for s in mundo["sitios"]}
    ok = {s: r for s, r in res.items() if r["costo"] is not None}
    mejor = None
    if ok:
        mejor_cob = max(r["cob"] for r in ok.values())
        cand = {s: r for s, r in ok.items() if r["cob"] >= mejor_cob - TOL_COB * n_total}
        mejor = min(cand, key=lambda s: cand[s]["costo"])
    if verbose:
        print(f"Semilla {semilla}: zonas con 100/500/1000 Mbps =", [list(mz.values()).count(v) for v in NIVELES])
        for s in mundo["sitios"]:
            r = res[s]
            print(f"  sitio {s}: cobertura {100 * r['cob'] / n_total:5.1f}% | estado {r['estado']:10s} |",
                  "costo ---" if r["costo"] is None else f"costo {r['costo']:8.1f} | puertos {r['puertos']:3d} | t {r['t']:.1f}s")
        print("  -> mejor sitio:", mejor)
    return mejor, ok.get(mejor)

if __name__ == "__main__":
    fuente = sys.argv[1] if len(sys.argv) > 1 else "sintetico"
    semillas = [int(a) for a in sys.argv[2:]] or [1]
    mundo = mundo_sintetico() if fuente == "sintetico" else cargar_mundo(fuente)
    for semilla in semillas:
        mejor, resultado = barrido(semilla, mundo)
        if mejor is None:
            continue
        ruta_json = f"solucion_pon_semilla{semilla}.json"
        with open(ruta_json, "w", encoding="utf-8") as archivo:
            json.dump(dict(semilla=semilla, sitio=mejor, resultado=resultado), archivo)
        ruta_png = f"red_pon_semilla{semilla}.png"
        graficar(mundo, resultado, mejor, semilla, ruta_png, mostrar=False)
        print("Guardado:", ruta_json, ruta_png)