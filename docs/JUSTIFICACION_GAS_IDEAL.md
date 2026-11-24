# Justificación del gas ideal (neón / helio)

Documento de apoyo para el PDF y la presentación. Explica por qué se puede
modelar la expansión con un gas ideal, por qué un gas noble monoatómico es una
elección razonable y cómo se mapean las unidades lattice a magnitudes físicas.

---

## 1. Objetivo y criterio

El objetivo del laboratorio es mostrar **paralelismo híbrido** con un autómata
celular, no reproducir con precisión la termodinámica del neón. Por eso se
busca el modelo más simple que:

- tenga reglas **locales** y un estado finito por celda,
- conserve **masa** y **momento**,
- produzca una expansión con morfología razonable,
- y sea defendible físicamente con supuestos explícitos.

El gas ideal cumple esos requisitos. Las correcciones de gas real se cuantifican
abajo y resultan despreciables a la presión del experimento.

---

## 2. Supuestos del modelo

1. **Partículas idénticas y puntuales**, sin estructura interna.
2. **Sin fuerzas de largo alcance** (nada de Lennard-Jones, Coulomb, etc.).
3. **Colisiones binarias elásticas y locales**, que conservan masa y momento.
4. **Caos molecular**: las velocidades de partículas que colisionan no están
   correlacionadas (hipótesis estándar para derivar las ecuaciones de
   transporte).
5. **Gas diluido o moderadamente diluido**: la distancia entre partículas es
   grande frente al alcance de las interacciones.
6. **Gas monoatómico**: sin grados de libertad rotacionales ni vibracionales;
   `γ = cp/cv = 5/3`.
7. **Aproximación isotérmica** en el autómata: la temperatura no es una
   variable dinámica (ver §5).

En el autómata, la “colisión” no es un potencial intermolecular: es una regla
local que redistribuye direcciones conservando los invariantes. La **exclusión**
(máximo una partícula por dirección por celda) es una restricción de ocupación
del lattice, no una fuerza física; sirve para mantener el estado finito y para
que las colisiones sean localmente binarias.

---

## 3. Por qué neón o helio

Los gases nobles monoatómicos son los candidatos más cercanos a un gas ideal:

- Capa electrónica cerrada: interacciones de van der Waals muy débiles.
- Átomo esférico: no hay momentos dipolares permanentes.
- Sin modos internos: `γ = 5/3`.
- A presión atmosférica y temperatura ambiente, su factor de compresibilidad es
  prácticamente 1.

| Propiedad | Neón | Helio |
|---|---:|---:|
| Masa molecular `m` | 20.1797 u = 3.3509×10⁻²⁶ kg | 4.0026 u = 6.6465×10⁻²⁷ kg |
| `γ` | 5/3 | 5/3 |
| `c_s` isotérmica a 300 K | 351.6 m/s | 789.4 m/s |
| `c_s` adiabática a 300 K | 453.9 m/s | 1019.1 m/s |
| Densidad a 1 atm, 300 K | 0.8197 kg/m³ | 0.1626 kg/m³ |
| Segundo coeficiente virial `B2` (300 K) | ≈ 11.4 cm³/mol | ≈ 11.8 cm³/mol |

Se elige **neón** como caso base por sus escalas: a 1 atm y 300 K la velocidad
del sonido es moderada (352 m/s isotérmica, 454 m/s adiabática) y la densidad
es fácil de verificar (≈ 0.82 kg/m³). El helio es igual de válido; solo cambia
la escala temporal (es ~2.2 veces más rápido) y la masa total.

---

## 4. ¿Qué tan bueno es el gas ideal a 1 atm?

La desviación se mide con el factor de compresibilidad:

```
Z = p / (n k_B T) = 1 + B2 n + B3 n² + ...
```

A `T = 300 K` y `p = 1 atm = 101325 Pa`:

```
n = p / (k_B T) = 2.446×10²⁵ moléculas/m³ = 40.62 mol/m³
```

Con `B2(Ne) ≈ 11.4 cm³/mol = 11.4×10⁻⁶ m³/mol`:

```
B2 n ≈ 4.6×10⁻⁴
Z ≈ 1.00046
```

Es decir, **la corrección es del orden del 0.05%**. A 10 atm sería ~0.5% y a
100 atm ~5%; el experimento trabaja a 1 atm o menos, así que la aproximación de
gas ideal es excelente. Para helio la conclusión es la misma (su `B2` es del
mismo orden).

La exclusión del lattice no cambia esta conclusión: a densidad media
`ρ = 3` partículas/celda (media ocupación de las 6 direcciones) las colisiones
siguen siendo locales y binarias; la exclusión impone un techo de densidad,
pero no introduce una fuerza atractiva ni un potencial de largo alcance.

---

## 5. Del gas ideal al autómata lattice

En el gas ideal:

```
p = n k_B T
c_s,iso = sqrt(k_B T / m)          (isotérmica)
c_s,ad  = sqrt(γ k_B T / m)        (adiabática)
```

En el LGA de 6 direcciones y velocidad unitaria:

```
p = ρ / 3
c_s² = 1/3
```

donde `ρ` es la densidad en partículas/celda y `p` en unidades lattice. Es un
modelo **isotérmico**: todas las partículas tienen la misma rapidez y no hay
ecuación de energía. En el límite hidrodinámico, las reglas de colisión
conservan masa y momento y reproducen las ecuaciones de Navier–Stokes con una
ecuación de estado de gas ideal (más correcciones propias del lattice).

Consecuencia: el modelo **no** captura el enfriamiento adiabático de una
expansión real. Si el profesor pidiera una expansión adiabática, habría que
usar un LGA térmico (partícula en reposo, niveles de energía) o un LBM con
ecuación de energía. Para este laboratorio, la aproximación isotérmica es
suficiente y se declara explícitamente.

---

## 6. Qué se pierde al despreciar las interacciones

- Coeficientes de transporte reales (viscosidad, conductividad térmica).
- Cambios de temperatura durante la expansión.
- Correcciones del virial a alta presión.
- Estructura fina de choques y efectos de gas real.

Nada de esto es el objeto del laboratorio. La morfología de la expansión, la
conservación de invariantes y el comportamiento del paralelismo son robustos
frente a estas simplificaciones.

---

## 7. Mapeo de unidades lattice → SI

Se fija:

- `T0 = 300 K`
- `p0 = 1 atm = 101325 Pa`
- `Δx` = tamaño de celda en metros (se elige)
- `ρ0 = 3` partículas/celda ↔ `p0`, `T0`

Entonces:

```
n0        = p0 / (k_B T0)                 [moléculas/m³]
ρ_mass,0  = n0 · m                        [kg/m³]
c_s,iso   = sqrt(k_B T0 / m)              [m/s]
Δt        = (1/√3) · Δx / c_s,iso         [s]
v_lattice = Δx / Δt = √3 · c_s,iso        [m/s]
```

Si se elige un `Δx` macroscópico, cada partícula lattice representa un
conglomerado de moléculas:

```
N_cluster = n0 · Δx³ / 3
m_particle = N_cluster · m
```

Si se exige que **una partícula lattice = una molécula**, se despeja:

```
Δx = (3 k_B T0 / p0)^(1/3) ≈ 4.97 nm
```

### Tabla de mapeos

**Mapeo macroscópico** (`Δx = 1 mm`, `N = 256`):

| Magnitud | Neón | Helio |
|---|---:|---:|
| `Δx` | 1 mm | 1 mm |
| Dominio | 0.256 m | 0.256 m |
| `Δt` | 1.642 µs | 0.731 µs |
| `v_lattice` | 608.9 m/s | 1367.3 m/s |
| Moléculas por partícula | 8.154×10¹⁵ | 8.154×10¹⁵ |
| Masa por partícula | 2.732×10⁻¹⁰ kg | 5.420×10⁻¹¹ kg |
| Masa total (esfera, 412365 partículas) | 1.127×10⁻⁴ kg (0.113 g) | 2.235×10⁻⁵ kg (0.022 g) |

**Mapeo atómico** (1 partícula = 1 molécula):

| Magnitud | Neón | Helio |
|---|---:|---:|
| `Δx` | 4.968 nm | 4.968 nm |
| Dominio (`N=256`) | 1.272 µm | 1.272 µm |
| `Δt` | 8.159 ps | 3.634 ps |
| Masa total (esfera) | 1.382×10⁻²⁰ kg | 2.741×10⁻²¹ kg |

En ambos mapeos la dinámica lattice es idéntica; solo cambia la interpretación
de las escalas.

---

## 8. Conversión de las salidas del código

El `.vti` exporta `density = ρ` en partículas/celda (0..6). Para obtener
magnitudes físicas en ParaView o en un postproceso:

```
p [Pa]        = (ρ / 3) · p0
ρ_mass [kg/m³] = (ρ / 3) · ρ_mass,0
```

Ejemplo con neón, `p0 = 101325 Pa`, `ρ_mass,0 = 0.8197 kg/m³`:

- `density = 3` → `p = 1 atm`, `ρ_mass = 0.8197 kg/m³`
- `density = 0.3` → `p = 0.1 atm`, `ρ_mass = 0.0820 kg/m³`

En ParaView: `Calculator` → `density/3*101325` para presión en Pa.

---

## 8.1 Cámara de dos presiones

Para tener un gradiente de presión explícito y una razón finita, el código
soporta `GEOM=chamber`:

```
ρ_H = 3.0   → p_H = (3/3)·p0 = 1 atm
ρ_L = 0.3   → p_L = (0.3/3)·p0 = 0.1 atm
P_RATIO = 10
interfaz en x = IFACE_FRAC·N (default 0.25N)
RAMP = ancho del gradiente en celdas (0 = escalón)
```

A `t=0` desaparece la membrana y el gas de alta presión se expande hacia la
cámara de baja presión. Con `RAMP=8` el gradiente inicial es
`(p_H - p_L)/(8·Δx) ≈ 0.1125 atm/mm` para `Δx=1 mm`.

La masa total de la cámara es mayor que la de la esfera porque la cámara de
baja presión también contiene gas: para `N=256` y `ρ_H=3`, `ρ_L=0.3` es
≈ 16.36 millones de partículas lattice (≈ 4.47 g de neón con `Δx=1 mm`).

---

## 9. Cómo se selecciona en el código

Variables de entorno (el código imprime el mapeo al inicio):

```bash
GAS=neon DX_MM=1 T0_K=300 P0_ATM=1 \
  OMP_NUM_THREADS=4 mpirun -np 4 --bind-to none ./build/gas 256 1000 4 -1
```

- `GAS=neon|helium` (default `neon`)
- `DX_MM` = tamaño de celda en mm (default `1`); `DX_MM=0` activa el mapeo
  atómico (1 partícula = 1 molécula).
- `T0_K` (default `300`)
- `P0_ATM` (default `1`)

Salida típica:

```
gas=neon m=3.35092e-26 kg gamma=1.66667
T0=300 K p0=101325 Pa dx=0.001 m
cs_iso=351.577 m/s cs_ad=453.884 m/s dt=1.64217e-06 s v_lat=608.949 m/s
rho_mass=0.81974 kg/m3 molecules_per_particle=8.15438e+15
domain=0.256 m mass_total=0.000112677 kg
```

---

## 10. Referencias

- Dymond, J. H. y Smith, E. B., *The Virial Coefficients of Pure Gases and
  Mixtures*, Oxford University Press.
- Chapman, S. y Cowling, T. G., *The Mathematical Theory of Non-Uniform Gases*,
  Cambridge University Press.
- Frisch, U., Hasslacher, B. y Pomeau, Y., “Lattice-Gas Automata for the
  Navier-Stokes Equation”, *Phys. Rev. Lett.* 56, 1505 (1986).
- CRC Handbook of Chemistry and Physics (propiedades del neón y el helio).
- Documentación del curso en `ejemplos/`.
