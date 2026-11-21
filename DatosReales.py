# === Librerías ===
import json, math, os
import numpy as np, networkx as nx, geopandas as gpd
from pyproj import Transformer
from scipy.spatial import cKDTree

# === Parámetros (ajustar; los marcados con * son supuestos de diseño) ===
COMUNAS = ["Concepción", "Hualpén", "San Pedro de la Paz", "Talcahuano", "Chiguayante", "Penco"]
CENTRO_DIR = "Plaza de la Independencia, Concepción, Chile"   # * centro para elegir sitios candidatos de OLT
N_TOTAL = 3000                    # * clientes totales (caso base)
DMAX_KM = 20.0                    # alcance máximo por calle
CELDA_M = 1500                    # * tamaño de celda de demanda (m): cada celda es un punto de demanda
ZONA_M = 4500                     # * tamaño de zona de distribución (m): agrupa celdas para el sorteo de tasas
J_M = 1500                        # * separación de nodos candidatos de splitter (m)
SITIOS_M, SITIOS_RADIO_KM = 3000, 6.0   # * grilla y radio (km) alrededor del centro para sitios candidatos de OLT
RADIO_DIST_KM = 3.0               # * distancia máxima splitter -> punto de demanda considerada
RIESGO = "riesgo.geojson"         # opcional: polígonos donde no se permite OLT ni splitters (capas oficiales de riesgo)
SALIDA = "datos_pon.json"

# === 1. Descarga desde OpenStreetMap (requiere internet) ===
def descargar():
    import osmnx as ox
    comunas = ox.geocode_to_gdf([f"{c}, Región del Biobío, Chile" for c in COMUNAS])
    zona = comunas.geometry.union_all().buffer(0.01)             # ~1 km de holgura
    G = ox.graph_from_polygon(zona, network_type="drive")
    edif = ox.features_from_polygon(zona, tags={"building": True})
    centro = ox.geocode(CENTRO_DIR)                                # (lat, lon)
    return G, edif, centro, comunas

# === 2. Procesamiento (no usa internet: se puede probar con datos simulados) ===
def procesar(G, edif, centro, comunas=None, riesgo=None):
    # grafo no dirigido simple con largo en metros
    H = nx.Graph()
    for u, v, d in G.edges(data=True):
        l = float(d["length"])
        if H.has_edge(u, v): l = min(l, H[u][v]["length"])
        H.add_edge(u, v, length=l)
    lon0 = float(np.mean([G.nodes[v]["x"] for v in G]))
    epsg = 32700 + int((lon0 + 180) // 6) + 1                     # UTM sur
    tr = Transformer.from_crs(4326, epsg, always_xy=True)
    nodos = list(H.nodes)
    xy = np.array([tr.transform(G.nodes[v]["x"], G.nodes[v]["y"]) for v in nodos])
    arbol = cKDTree(xy)
    nodo_mas_cercano = lambda x, y: nodos[int(arbol.query([x, y])[1])]

    # zonas de riesgo (opcional): se excluyen nodos dentro de los polígonos
    prohibidos = set()
    if riesgo is not None and len(riesgo):
        from shapely.geometry import Point
        union = riesgo.to_crs(epsg).geometry.union_all()
        prohibidos = {v for v, (x, y) in zip(nodos, xy) if union.contains(Point(x, y))}

    # demanda: edificios -> celdas; clientes proporcionales al área construida (aproximación)
    e = edif[edif.geometry.type.isin(["Polygon", "MultiPolygon"])].to_crs(epsg)
    cen = e.geometry.centroid
    area = e.geometry.area.values
    ix, iy = np.floor(cen.x.values / CELDA_M).astype(int), np.floor(cen.y.values / CELDA_M).astype(int)
    celdas = {}
    for a_, cx, cy, p in zip(area, cen.x.values, cen.y.values, zip(ix, iy)):
        c = celdas.setdefault(p, [0.0, 0.0, 0.0]); c[0] += a_; c[1] += a_ * cx; c[2] += a_ * cy
    tot = sum(c[0] for c in celdas.values())
    cruda = {p: N_TOTAL * c[0] / tot for p, c in celdas.items()}
    n = {p: int(math.floor(v)) for p, v in cruda.items()}
    for p in sorted(cruda, key=lambda p: cruda[p] - n[p], reverse=True)[: N_TOTAL - sum(n.values())]: n[p] += 1
    k = max(1, ZONA_M // CELDA_M)
    puntos, nodo_punto = {}, {}
    for p, q in n.items():
        if q <= 0: continue
        pid = f"p{p[0]}_{p[1]}"
        c = celdas[p]
        puntos[pid] = dict(n=q, zona=f"z{p[0] // k}_{p[1] // k}")
        nodo_punto[pid] = nodo_mas_cercano(c[1] / c[0], c[2] / c[0])

    # nodos candidatos de splitter (uno por celda de J_M) y sitios de OLT
    def uno_por_celda(paso, filtro=lambda x, y: True):
        mejor = {}
        for v, (x, y) in zip(nodos, xy):
            if v in prohibidos or not filtro(x, y): continue
            c = (int(x // paso), int(y // paso))
            d = math.hypot(x - (c[0] + 0.5) * paso, y - (c[1] + 0.5) * paso)
            if c not in mejor or d < mejor[c][0]: mejor[c] = (d, v)
        return [v for _, v in mejor.values()]
    cx_, cy_ = tr.transform(centro[1], centro[0])
    J = uno_por_celda(J_M)
    sitios = uno_por_celda(SITIOS_M, lambda x, y: math.hypot(x - cx_, y - cy_) <= SITIOS_RADIO_KM * 1000)

    # distancias por calle
    nm, ek = str, (lambda a, b: "|".join(sorted((str(a), str(b)))))
    largo = {ek(a, b): d["length"] / 1000 for a, b, d in H.edges(data=True)}
    por_nodo = {}
    for pid, v in nodo_punto.items(): por_nodo.setdefault(v, []).append(pid)
    dsj, rutas, dji = {}, {}, {}
    for s in sitios:
        dist, paths = nx.single_source_dijkstra(H, s, cutoff=DMAX_KM * 1000, weight="length")
        dsj[nm(s)] = {nm(j): dist[j] / 1000 for j in J if j in dist}
        rutas[nm(s)] = {nm(j): [ek(a, b) for a, b in zip(paths[j], paths[j][1:])] for j in J if j in dist}
    for j in J:
        dd = nx.single_source_dijkstra_path_length(H, j, cutoff=RADIO_DIST_KM * 1000, weight="length")
        dji[nm(j)] = {pid: dd[v] / 1000 for v, pids in por_nodo.items() if v in dd for pid in pids}
    mundo = dict(n={p: d["n"] for p, d in puntos.items()}, zona={p: d["zona"] for p, d in puntos.items()},
                 J=[nm(j) for j in J], sitios=[nm(s) for s in sitios], dsj=dsj, dji=dji, rutas=rutas, largo=largo)

    # verificación del límite de 20 km: distancia por calle de cada sitio al centro de cada comuna
    if comunas is not None:
        print("Distancia por calle desde cada sitio al centro de cada comuna (km); '-' = fuera de 20 km")
        cc = comunas.to_crs(epsg).geometry.centroid
        nodos_com = [nodo_mas_cercano(p.x, p.y) for p in cc]
        print("sitio".ljust(14), *[str(c)[:10].ljust(11) for c in comunas.get("name", comunas.index)])
        for s in sitios:
            dist = nx.single_source_dijkstra_path_length(H, s, cutoff=DMAX_KM * 1000, weight="length")
            print(nm(s).ljust(14), *[(f"{dist[v] / 1000:.1f}" if v in dist else "-").ljust(11) for v in nodos_com])
    print(f"puntos de demanda: {len(puntos)} | zonas: {len({d['zona'] for d in puntos.values()})} | "
          f"nodos J: {len(J)} | sitios: {len(sitios)} | clientes: {sum(mundo['n'].values())}")
    return mundo

if __name__ == "__main__":
    G, edif, centro, comunas = descargar()
    riesgo = gpd.read_file(RIESGO) if os.path.exists(RIESGO) else None
    mundo = procesar(G, edif, centro, comunas, riesgo)
    with open(SALIDA, "w", encoding="utf-8") as f: json.dump(mundo, f)
    print("Guardado:", SALIDA)