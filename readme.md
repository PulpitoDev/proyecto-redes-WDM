
# Diseño y optimización de una red PON

## 1. Alcance del proyecto

El proyecto aborda el diseño y la optimización de una red de acceso PON (Passive Optical Network), considerando la selección de una ubicación óptima para la OLT y la distribución de splitters para atender la demanda de clientes.

El alcance se limita exclusivamente a la red PON. Las redes troncales quedan fuera del análisis.

### 1.1. Área geográfica

El estudio considera las siguientes comunas de la Región del Biobío, Chile:

- Concepción
- Hualpén
- San Pedro de la Paz
- Talcahuano
- Chiguayante
- Penco

### 1.2. Fuente de datos

Se utilizan datos geográficos reales de OpenStreetMap (OSM):

- **Red vial:** utilizada para calcular las distancias de recorrido entre los distintos elementos de la red.
- **Edificios:** utilizados para estimar y distribuir espacialmente la demanda de clientes.

### 1.3. Dimensión del caso de estudio

| Elemento | Cantidad |
|---|---:|
| Clientes | 3.000 |
| Puntos de demanda | 127 |
| Tamaño de celda | 1,5 km |
| Zonas de servicio | 31 |
| Nodos candidatos a splitter | 240 |
| Sitios candidatos a OLT | 20 |

---

## 2. Función objetivo del proyecto

El objetivo general es seleccionar la ubicación de la OLT que permita maximizar la cobertura de clientes y, posteriormente, minimizar el costo total de implementación de la red PON.

Para ello, se evalúa cada uno de los 20 sitios candidatos a OLT mediante la resolución independiente del modelo de optimización.

### 2.1. Criterio de selección de la OLT

La selección se realiza mediante un criterio jerárquico:

1. **Maximización de cobertura:** se identifica la máxima cobertura alcanzable entre todos los sitios candidatos. Se consideran aceptables aquellos sitios cuya cobertura esté a menos de 0,5 puntos de la cobertura máxima.
2. **Minimización de costos:** entre los sitios aceptables, se selecciona aquel que presenta el menor costo total según la función objetivo `F(x)`.

### 2.2. Función objetivo

Para una ubicación determinada de la OLT `s`, el modelo busca minimizar:

```math
\begin{aligned}
\min F(x) ={}&
\sum_{j,c}
(q_t+k_r+\beta L(s,j))x(j,c)\\
&+\sum_e \alpha \ell(e)g(e)\\
&+\sum_{i,j,c}
\beta n(i)\delta(i,j)a(i,j,c)
\end{aligned}
```

La función objetivo considera tres componentes:

**1. Costos de equipamiento y fibra de alimentación**

```math
\sum_{j,c}(q_t+k_r+\beta L(s,j))x(j,c)
```

Considera el costo de los puertos, los splitters y la fibra de alimentación entre la OLT y los splitters seleccionados.

**2. Costos de zanja de alimentación**

```math
\sum_e \alpha \ell(e)g(e)
```

Representa el costo de construcción de zanjas para la infraestructura de alimentación. Cuando varios splitters comparten un mismo tramo, este se contabiliza una sola vez.

**3. Costos de fibra de distribución**

```math
\sum_{i,j,c}\beta n(i)\delta(i,j)a(i,j,c)
```

Representa el costo de la fibra de distribución entre los nodos candidatos a splitter y los puntos de demanda atendidos.

---

## 3. Restricciones del modelo

### 3.1. Asignación de demanda

Cada punto de demanda debe ser atendido por un único splitter.

### 3.2. Capacidad de los splitters

La cantidad de clientes asignados a un splitter no puede superar su razón de división:

```math
\text{Clientes del splitter} \leq r
```

Donde:

```math
r \in \{8,16\}
```

### 3.3. Requerimientos de tasa de transmisión

La tasa disponible por cliente debe satisfacer la tasa mínima exigida por su zona, tanto en bajada como en subida.

Se considera:

```math
\text{Tasa por cliente} =
\frac{\text{Tasa de línea}}{r}
```

Por lo tanto:

```math
\frac{T_{\text{línea}}}{r}
\geq T_{\min,z}
```

Donde `T_min,z` corresponde al requerimiento mínimo de transmisión de la zona `z`.

### 3.4. Restricción de alcance

La distancia total de la ruta entre la OLT, el splitter y el punto de demanda debe respetar el alcance máximo permitido:

```math
L(s,j)+\delta(i,j)
\leq \min(20\text{ km},L_{\max}(r))
```

Se consideran los siguientes límites:

| Razón de división | Alcance máximo |
|---|---:|
| 1:8 | 20 km |
| 1:16 | 14 km |

### 3.5. Capacidad de la OLT

La cantidad total de puertos utilizados en la OLT no puede superar los 400 puertos disponibles.

```math
\text{Puertos totales} \leq 400
```

### 3.6. Compartición de infraestructura

Los tramos de zanja compartidos por diferentes splitters se contabilizan una única vez, evitando duplicar los costos de infraestructura común.

---

## 4. Supuestos del modelo

### 4.1. Supuestos generales del caso

| Parámetro | Supuesto |
|---|---|
| Número de clientes | 3.000 |
| Cantidad de OLT | Una única OLT |
| Tecnología | PON |
| Respaldo de red | No considerado |
| Zonas de riesgo | No consideradas |
| Enlace de salida de la OLT | Se asume existente |

### 4.2. Supuestos de demanda

- Los clientes se distribuyen proporcionalmente al área construida de los edificios identificados en OpenStreetMap, y no según datos censales.
- La demanda total corresponde a 3.000 clientes.
- Cada zona tiene un requerimiento mínimo de tasa de transmisión.
- La tasa mínima por zona se selecciona aleatoriamente entre 100, 500 y 1.000 Mbps, con igual probabilidad.
- Los requerimientos de bajada y subida son simétricos.

### 4.3. Supuestos de red

- La tasa efectiva por cliente se calcula como la tasa de línea dividida por la razón de división del splitter.
- No se considera overhead de transmisión.
- Se supone que todos los clientes están activos simultáneamente.
- Se permiten razones de división de 1:8 y 1:16.
- Se considera una única etapa de splitters.
- Se asigna una fibra de alimentación por splitter.
- Se considera una fibra de distribución por cliente.

### 4.4. Parámetros de costos y capacidades

Los siguientes parámetros son valores ficticios utilizados para la experimentación y evaluación del modelo. No corresponden a cotizaciones comerciales ni a fuentes oficiales.

| Parámetro | Valores considerados |
|---|---|
| Costo de puerto OLT | 10 / 18 / 25 |
| Costo de splitter | 0,8 / 1,2 |
| Costo de fibra | 0,5 por km |
| Costo de zanja | 8,0 por km |
| Alcance máximo para 1:8 | 20 km |
| Alcance máximo para 1:16 | 14 km |
| Capacidad máxima de puertos OLT | 400 |

### 4.5. Supuestos de discretización espacial

Para reducir la complejidad computacional del problema, se utiliza una representación espacial discretizada:

| Parámetro | Valor |
|---|---:|
| Tamaño de las celdas de demanda | 1,5 km |
| Separación entre nodos candidatos | 1,5 km |
| Distancia máxima splitter → punto de demanda | 3 km |

---

## 5. Limitaciones del modelo

Los resultados obtenidos deben interpretarse considerando las siguientes limitaciones:

- **Estimación de demanda:** los clientes se distribuyen en función del área construida y no de información censal o registros reales de abonados.
- **Discretización espacial:** la demanda se agrupa en celdas de 1,5 km, lo que puede introducir diferencias entre la ubicación modelada y la distribución real de los clientes.
- **Ubicación de splitters:** las posiciones de los nodos candidatos son aproximadas, con una precisión espacial del orden de 1 km.
- **Infraestructura existente:** no se modela la infraestructura de red previamente desplegada.
- **Respaldo:** no se consideran rutas alternativas ni mecanismos de redundancia.
- **Riesgo:** no se incorporan restricciones asociadas a zonas de riesgo o vulnerabilidad de infraestructura.
- **Distancias:** los recorridos se calculan utilizando la red vial de OpenStreetMap, sin considerar pasajes peatonales.
- **Alcance:** la restricción de 20 km se evalúa hasta el nodo representativo de la celda de demanda.
- **Parámetros económicos:** los costos utilizados son ficticios, por lo que los resultados dependen directamente de los valores definidos para puertos, splitters, fibra y zanjas.

---

## 6. Consideraciones finales

El modelo permite abordar el diseño de una red PON desde una perspectiva de optimización geográfica y económica, integrando información territorial real de OpenStreetMap con restricciones técnicas de capacidad, cobertura, tasa de transmisión y alcance.

La metodología permite comparar diferentes ubicaciones candidatas para la OLT, priorizando la cobertura y minimizando los costos de despliegue de la red de acceso.

Es importante destacar que los resultados corresponden a una evaluación de diseño bajo los supuestos definidos, y no necesariamente a un presupuesto real de implementación ni a un diseño constructivo definitivo.
