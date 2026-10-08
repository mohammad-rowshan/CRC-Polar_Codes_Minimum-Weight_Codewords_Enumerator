# Minimum-weight codewords of polar, PAC, PS-PAC and CRC-polar codes

This directory contains two Python implementations of the same calculation:

- `mwc_enumerator.py` uses the Python standard library only.
- `mwc_enumerator_numba.py` uses Numba to accelerate the enumeration loops.

Both scripts consider **one code at a time**. The code family is selected by
`code_type`: `polar`, `pac`, `ps-pac`, or `crc-polar`. The code parameters are
specified by the user; the minimum distance and its multiplicity are calculated
by the script. No target weight is required.

The usual calculation uses the structure of minimum-weight codewords of polar
codes and its extension to pre-transformed polar codes [1, 3-5]. When this
calculation does not establish the minimum distance, the script can search for
low-weight codewords by means of parity-check syndromes. That additional search
is described in Section 5; it is not part of every run.

## 1. Running the scripts

Python 3.10 or later is required. The standard Python script has no additional
dependencies. For the accelerated script, install Numba and NumPy:

```bash
python -m pip install numpy numba
```

The usual way to run a code is to edit `CONFIG` near the beginning of either
script. For example, the configuration for a PS-PAC (64,32) code is

```python
CONFIG = dict(
    code_type="ps-pac",
    N=64,
    K=32,
    construction="dega",
    design_snr_db=2.0,
    alpha=7,
    pac_poly="1011011111",
)
```

Save the file and run either version:

```bash
python mwc_enumerator.py
```

or

```bash
python mwc_enumerator_numba.py
```

The scripts can also be run from the command line. For example,

```bash
python mwc_enumerator.py --code-type ps-pac --N 64 --K 32 --construction dega --design-snr-db 2 --alpha 7 --pac-poly 1011011111
```

Use `mwc_enumerator_numba.py` in place of `mwc_enumerator.py` to run the
accelerated version. **Command-line parameters replace `CONFIG` entirely**;
unspecified options then take their command-line defaults. The available
options can be listed with `python mwc_enumerator.py --help`.

The function `enumerate_code(...)` may also be called from another Python
program. It returns the same results shown on the screen.

## 2. Code parameters

`N` denotes the block length and `K` denotes the number of independent
information bits. For CRC-polar codes, `K` excludes the CRC bits: a CRC of
degree `t` gives `K + t` selected polar input positions.

The information set is chosen with `construction`:

| Value | Construction |
| --- | --- |
| `dega` | Gaussian approximation (GA), with `design_snr_db` and optionally `design_rate` [7] |
| `crossed-rm` | Crossed RM-polar rate profile of [6], using GA reliability within the boundary row-weight class |
| `5g` | Reliability sequence from 3GPP TS 38.212 [8] |
| `pw` | Polarization-weight construction; `pw_beta` can be adjusted [9, 10] |
| `explicit` | Information indices supplied through `info_set` |

For example, `info_set="15,23,26-31,38,39"` specifies an explicit information
set. The 5G option selects positions from the standard reliability sequence
for `N <= 1024`; it does not perform rate matching or interleaving.

All indices begin at zero. The polar transform is

$$
G_N = F^{\otimes n}, \qquad
F = \begin{bmatrix}1&0\\1&1\end{bmatrix}, \qquad N=2^n.
$$

No separate bit-reversal permutation is applied.

The other main parameters are:

| Parameter | Meaning |
| --- | --- |
| `pac_poly` | Convolution coefficients in increasing delay order, e.g. `1011011` represents `(1,0,1,1,0,1,1)` |
| `alpha` | Number of exchanged positions in PS-PAC |
| `crc_poly` | CRC polynomial, e.g. `x^5+x^3+1` or `0x29` |
| `crc_positions` | CRC positions; by default the last `t` selected positions |
| `per_coset` | Print the contribution of each leading coset |

The PAC precoder is defined by

$$
u_i = \sum_{\ell=0}^{s} p_{\ell} v_{i-\ell} \pmod 2,
\qquad v_j=0 \text{ for } j<0.
$$

Here `v` is the rate-profile vector and `u` is the input to the polar
transform. For CRC-polar codes, the remainder is computed from

$$
r(x) = d(x)x^t \bmod q(x).
$$

CRC positions must follow the data positions in the selected information set.
For PS-PAC with `construction="explicit"`, `info_set` is understood to be
**already shifted**, unless `reliability_order` is also supplied. With an
explicit reliability order, the shift is applied to the specified base set.

## 3. Minimum-weight enumeration

### 3.1. Polar codes

The polar calculation uses the formula of Rowshan, Dau and Viterbo [1]. For a
polar information set satisfying the standard partial-order property, the
weight of row `i` of the polar transform is

$$
w(g_i) = 2^{\mathrm{wt}(i)},
$$

where `wt(i)` is the Hamming weight of the binary expansion of `i`. Define

$$
\begin{aligned}
d &= \min_{i\in\mathcal I} w(g_i),\\
\mathcal B &= \lbrace i\in\mathcal I : w(g_i)=d \rbrace,\\
\mathcal K_i &= \lbrace j\in\mathcal I : j>i,\ |
    \mathrm{supp}(j)\setminus\mathrm{supp}(i)|=1 \rbrace.
\end{aligned}
$$

Then

$$
d_{\min}=d, \qquad
A_{d_{\min}} = \sum_{i\in\mathcal B} 2^{|\mathcal K_i|}.
$$

The script constructs these sets and evaluates the sum. The minimum-weight
count is equivalent to the earlier result for decreasing monomial codes by
Bardet, Dragoi, Otmani and Tillich [2]. The implementation follows [1]. It
does not use the formula as an exact count when the required partial-order
condition is not satisfied.

### 3.2. PAC codes

The PAC method follows the formation of polar minimum-weight codewords [1],
the fast enumeration method in [3], and the capable-coset analysis in [4].

A polar minimum-weight word with leading row `i` is formed using indices

$$
\mathcal O = \lbrace i\rbrace\cup\mathcal J\cup\mathcal M(\mathcal J),
\qquad \mathcal J\subseteq\mathcal K_i,
$$

where the balancing set `M(J)` ensures that the weight remains unchanged.
When a new element is included in `J`, the balancing set is updated by the
`addToM` rule of [3, 4]; repeated terms cancel over the binary field.

For each leading coset, frozen indices that can reduce its minimum-weight
multiplicity are identified by

$$
\mathcal F_i^* = \lbrace f\in\mathcal I^c : f>i,\ |
\mathrm{supp}(f)\setminus\mathrm{supp}(i)|>1 \rbrace.
$$

When this set is empty, the coset is **incapable** of reducing the polar
minimum-weight count [4]. Otherwise, let `f_max` be its largest element. The
script enumerates subsets of the truncated core

$$
\mathcal K_i^{\mathrm{tr}} = \lbrace j\in\mathcal K_i : j<f_{\max} \rbrace.
$$

For a given subset, the required polar input bits are compared with the bits
obtained by convolutional precoding. At frozen positions `v_i=0`, so any
inconsistent subset is discarded. The choices beyond the truncated core are
accounted for by a multiplicity factor, as in [3, 4]. The surviving counts
are added over the leading cosets.

Zunker, Geiselhart and ten Brink [13] proposed a tree-intersection method
for enumerating minimum-weight codewords of pre-transformed polar codes.
Their method is not implemented here.

This is not an enumeration of all `2^K` information vectors. The parameter
`max_structural_subsets` (default `16777216`) limits the number of subsets
examined in an individual truncated-core calculation. This limit concerns the
structural PAC calculation, not the syndrome search in Section 5.

The method counts a specified polar-row weight class. When that information
alone does not determine the minimum distance of the PAC code, further
low-weight checking is required.

### 3.3. Profile-shifted PAC (PS-PAC) codes

PS-PAC codes were introduced by Gu, Rowshan and Yuan [5]. The information set
is modified before convolutional precoding: the `alpha` most reliable
positions of the base information set are frozen, while an equal number of
positions are admitted from the frozen set according to the reliability
ordering. The shifted information set, denoted `I'`, still has `K` elements.

The convolutional precoder is applied to `I'`, not to the original set. For a
newly frozen position `f`, the profile bit satisfies `v_f=0`, but the polar
input bit may be nonzero:

$$
u_f = \sum_{\ell=1}^{s} p_{\ell}v_{f-\ell}\pmod 2.
$$

Such constraint-dependent positions can prevent the formation of
minimum-weight codewords [5]. The script therefore constructs the shifted
profile first, then applies the PAC coset enumeration with the new frozen
positions. It does **not** calculate an ordinary PAC multiplicity and adjust
that number afterward.

A profile shift can introduce generator rows of smaller weight. Even if their
structural minimum-weight contributions vanish, this alone does not exclude
codewords of intermediate weights. Section 5 explains the additional check
used in that case. Ellouze *et al.* [12] studied the minimum-distance
properties of PAC codes using polar cosets and prefix pruning. Their
algorithm is not implemented here.

### 3.4. CRC-polar codes

For CRC-polar codes, the script uses the polar minimum-weight construction
[1] together with the CRC constraints studied by Gu, Rowshan and Yuan [5].
Let `I_CRC` consist of the `K` data positions and `t` CRC positions `R`. The
minimum generator-row weight of the underlying polar code is

$$
w_0 = \min_{i\in\mathcal I_{\mathrm{CRC}}} w(g_i).
$$

For each possible free leading row of weight `w_0`, the script visits subsets
`J` of its set `K_i` and constructs the corresponding balancing set. The
candidate polar input has support

$$
\mathcal O = \lbrace i\rbrace\cup\mathcal J\cup\mathcal M(\mathcal J).
$$

A candidate is rejected if it uses an ordinary frozen position or if its
bits at the CRC positions disagree with the remainder calculated from the
data bits. Otherwise, the word is counted. Contributions from different
leading cosets are then added to obtain `A_w0` for the CRC-polar code.

The script uses direct subset traversal and CRC testing. It does not
implement every early-pruning step of Algorithm 1 in [5]. If `A_w0 > 0`,
then `w0` is the minimum distance. If `A_w0 = 0`, the CRC has eliminated that
weight class, but words of weight `1.5 w0` or other intermediate weights may
remain [5]. The separate calculation in Section 5 is then needed.

## 4. Examples

For each example, replace `CONFIG` in either script and run it without
command-line arguments. The two versions use the same configuration names.

### Polar (128,64), crossed RM-polar

```python
CONFIG = dict(
    code_type="polar", N=128, K=64,
    construction="crossed-rm", design_snr_db=3.5,
)
```

Expected: `d_min = 16`, `A_dmin = 94488`.

### PAC (128,64), crossed RM-polar

```python
CONFIG = dict(
    code_type="pac", N=128, K=64,
    construction="crossed-rm", design_snr_db=3.5,
    pac_poly="1011011",
)
```

Expected: `d_min = 16`, `A_dmin = 3120` [3, 4].

### PS-PAC (64,32), GA construction

```python
CONFIG = dict(
    code_type="ps-pac", N=64, K=32,
    construction="dega", design_snr_db=2.0,
    alpha=7, pac_poly="1011011111",
)
```

Expected: `d_min = 8`, `A_dmin = 23`, with no words of weights 4, 5, 6 or 7.
The low-weight check is required for this example. The polynomial above is
the one reported with `A_8 = 23` in the numerical table of [5]. Another
passage of the manuscript gives `1011011011`, which does not reproduce the
tabulated result.

### CRC-polar (64,32+5)

```python
CONFIG = dict(
    code_type="crc-polar", N=64, K=32,
    construction="dega", design_snr_db=2.0,
    crc_poly="x^5+x^3+1",
)
```

Expected: `d_min = 8`, `A_dmin = 64` [5].

### CRC-polar (128,64+8), polarization-weight profile

```python
CONFIG = dict(
    code_type="crc-polar", N=128, K=64,
    construction="pw", crc_poly="0x1D5",
)
```

Expected: `d_min = 8`, `A_dmin = 14`. This result is specific to the
information set and CRC polynomial shown.

### CRC-polar (32,16+5), increased minimum distance

```python
CONFIG = dict(
    code_type="crc-polar", N=32, K=16,
    construction="explicit", info_set="7,11-15,17-31",
    crc_poly="x^5+x^3+1",
)
```

The CRC removes all weight-four words. The additional low-weight search
finds `d_min = 6` and `A_dmin = 14`. Thus the minimum distance of a CRC-polar
code need not be a power of two.

## 5. Low-weight search using parity-check syndromes

**This search is used only when the structural calculation has not
established the true minimum distance.** It is not needed for the usual polar
formula or for CRC-polar codes whose lowest polar weight class survives.

The script obtains a generator matrix for the selected code, including the
convolutional precoder or CRC constraints, and derives a parity-check matrix
`H`. Denote its `j`-th column by `h_j`. A support `S` corresponds to a
codeword precisely when

$$
\bigoplus_{j\in S}h_j=0.
$$

Divide the coordinates into disjoint sets `L` and `R`. For any subset `T`,
write its syndrome as

$$
\sigma(T)=\bigoplus_{j\in T}h_j.
$$

Let `M_L(a,s)` be the number of subsets of `L` of size `a` and syndrome `s`;
define `M_R(b,s)` in the same way for `R`. A support with `a` positions in
`L` and `w-a` positions in `R` is a codeword if the syndromes of its two
parts agree. Consequently,

$$
A_w = \sum_{a=\max(0,w-|R|)}^{\min(w,|L|)}
      \sum_s M_L(a,s)M_R(w-a,s).
$$

The implementation forms a table of syndrome multiplicities for one part
and counts matches in the other. Distinct subsets may have the same
syndrome; their multiplicities must therefore be retained. Every support
has a unique split between `L` and `R`, so it is counted exactly once.
Testing successive weights determines `d_min` and `A_dmin`, provided all
smaller weights have been excluded.

For example, the even-parity `[3,2,2]` code has `H = [1 1 1]`. Choose
`L = {0}` and `R = {1,2}`. At weight two, the splits `(0,2)` and `(1,1)`
contribute one and two supports, respectively. Hence `A_2 = 3`.

This is a support enumeration rather than a search through `2^K` messages.
It uses standard parity-check algebra and a syndrome-collision procedure.
Stern's classical low-weight search [11] is related to the collision idea,
but Stern's randomized algorithm is **not** the exact weight-by-weight
counter implemented here. The expression above specifies the method used by
these scripts.

### Limits of the syndrome search

At large block lengths, the number of subsets can still be substantial.
The following parameters apply **only to this low-weight search**:

| Parameter | Default | Meaning |
| --- | ---: | --- |
| `max_exact_weight` | `10` | Highest weight examined by the syndrome search |
| `max_syndrome_combos` | `12000000` | Maximum subset combinations allowed for either side of a split |
| `max_hash_entries` | `1000000` | Maximum entries allowed in a syndrome table |

These are computational limits, **not user-supplied target weights**. They
may be changed in `CONFIG` or on the command line. If the search cannot
exclude all lower weights within the limits, the script reports
`NOT CERTIFIED` instead of assigning an unsupported value to `A_dmin`.
It may also return bounds on the minimum distance.

The Python result contains `d_min`, `A_dmin`, and `certified`. The
`weight_counts` field contains weights counted exactly by the syndrome
search. The `structural_classes` field contains counts obtained from
particular polar-row weight classes; those counts should not be interpreted
as a complete weight enumerator. Setting `per_coset=True` prints the
individual structural coset counts.

## 6. Execution time

Numba is useful when many subsets must be examined. For PAC (128,64) with
the crossed RM-polar profile and `pac_poly="1011011"`, an indicative
measurement on an x86-64 Linux system was approximately **1.9 s** with
pure Python and **0.13 s** with Numba after compilation. These are not
Windows measurements, and timings will depend on the processor and code
parameters.

**The first Numba run may take longer.** Numba uses just-in-time (JIT)
compilation: the relevant Python functions are converted to machine code
when they are first called. Subsequent calls normally use the compiled
functions. Caching is enabled where supported, although source changes
or different argument types can cause recompilation. The quoted Numba
time excludes initial compilation.

## References

[1] M. Rowshan, S. H. Dau and E. Viterbo,
"On the Formation of Min-Weight Codewords of Polar/PAC Codes and Its
Applications," *IEEE Transactions on Information Theory*, vol. 69,
no. 12, pp. 7627-7649, 2023.
[doi:10.1109/TIT.2023.3319015](https://doi.org/10.1109/TIT.2023.3319015).

[2] M. Bardet, V.-F. Dragoi, A. Otmani and J.-P. Tillich,
"Algebraic Properties of Polar Codes From a New Polynomial Formalism,"
*IEEE International Symposium on Information Theory (ISIT)*,
pp. 230-234, 2016.
[doi:10.1109/ISIT.2016.7541295](https://doi.org/10.1109/ISIT.2016.7541295).

[3] M. Rowshan and J. Yuan,
"Fast Enumeration of Minimum Weight Codewords of PAC Codes,"
*IEEE Information Theory Workshop (ITW)*, pp. 255-260, 2022.
[doi:10.1109/ITW54588.2022.9965901](https://doi.org/10.1109/ITW54588.2022.9965901).

[4] M. Rowshan and J. Yuan,
"On the Minimum Weight Codewords of PAC Codes: The Impact of
Pre-Transformation," *IEEE Journal on Selected Areas in Information
Theory*, vol. 4, pp. 487-498, 2023.
[doi:10.1109/JSAIT.2023.3312678](https://doi.org/10.1109/JSAIT.2023.3312678).

[5] X. Gu, M. Rowshan and J. Yuan,
"CRC-Polar-Inspired PAC Codes: Enumeration and Minimum-Weight Codewords
Pruning," accepted for publication in *IEEE Transactions on Communications*.

[6] M. Rowshan, A. Burg and E. Viterbo,
"Polarization-Adjusted Convolutional (PAC) Codes: Sequential Decoding vs
List Decoding," *IEEE Transactions on Vehicular Technology*, vol. 70,
no. 2, pp. 1434-1447, 2021.
[doi:10.1109/TVT.2021.3052550](https://doi.org/10.1109/TVT.2021.3052550).

[7] P. Trifonov, "Efficient Design and Decoding of Polar Codes,"
*IEEE Transactions on Communications*, vol. 60, no. 11,
pp. 3221-3227, 2012.
[doi:10.1109/TCOMM.2012.081512.110872](https://doi.org/10.1109/TCOMM.2012.081512.110872).

[8] 3GPP TS 38.212, *Multiplexing and Channel Coding*,
polar reliability sequence.

[9] Y. Zhou, R. Li, H. Zhang, H. Luo and J. Wang,
"Polarization Weight Family Methods for Polar Code Construction,"
*IEEE Vehicular Technology Conference (VTC Spring)*, 2018.
[doi:10.1109/VTCSpring.2018.8417498](https://doi.org/10.1109/VTCSpring.2018.8417498).

[10] G. He *et al.*,
"beta-expansion: A Theoretical Framework for Fast and Recursive Construction
of Polar Codes," arXiv:1704.05709, 2017.
[arXiv:1704.05709](https://arxiv.org/abs/1704.05709).

[11] J. Stern, "A Method for Finding Codewords of Small Weight,"
in *Coding Theory and Applications*, Lecture Notes in Computer Science,
vol. 388, pp. 106-113, Springer, 1989.
[doi:10.1007/BFb0019850](https://doi.org/10.1007/BFb0019850).

[12] M. Ellouze, R. Tajan, C. Leroux, C. Jego and C. Poulliat,
"Low-Complexity Algorithm for the Minimum Distance Properties of PAC Codes,"
*International Symposium on Topics in Coding (ISTC)*, 2023.
[doi:10.1109/ISTC57237.2023.10273572](https://doi.org/10.1109/ISTC57237.2023.10273572).

[13] A. Zunker, M. Geiselhart and S. ten Brink,
"Enumeration of Minimum Weight Codewords of Pre-Transformed Polar
Codes by Tree Intersection," *58th Annual Conference on Information
Sciences and Systems (CISS)*, pp. 1-6, 2024.
[doi:10.1109/CISS59072.2024.10480163](https://doi.org/10.1109/CISS59072.2024.10480163).
