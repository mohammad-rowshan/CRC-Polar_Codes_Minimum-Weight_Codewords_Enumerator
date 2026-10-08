"""Unified minimum-distance / MWC enumerator for polar, PAC, PS-PAC and CRC-polar.

ONE self-contained research script with embedded 3GPP 5G sequence.
Requires only Python >= 3.10 (no third-party packages). One code_type per run.
Edit CONFIG and run without arguments, or use command-line arguments.

Method: polar closed-form, structural PAC / CRC cosets for the leading weight,
then optional exact weight-limited syndrome search when that weight vanishes.
An incomplete resource-limited search NEVER reports a guessed A_dmin.

References: Rowshan--Dau--Viterbo (IEEE TIT 2023; equivalent polar
minimum-weight count to Bardet--Dragoi--Otmani--Tillich, ISIT 2016); Rowshan--Yuan
(ITW 2022, IEEE JSAIT 2023); Gu--Rowshan--Yuan (IEEE TCOM, accepted).
"""
from __future__ import annotations
import argparse
import itertools
import sys
import math
import re
import time
from dataclasses import dataclass
from typing import Iterable, Sequence

# --------------------- USER CONFIGURATION ------------------------------
# Used when the script is run without command-line arguments. You may also
# set CONFIG = None to require CLI arguments. Edit any keys as needed.
# Command-line arguments always override this configuration in its entirety.
CONFIG = dict(
    code_type="pac", N=128, K=64,
    construction="crossed-rm", design_snr_db=3.5,
    pac_poly="1011011",
)
"""CONFIG = dict(
    code_type="crc-polar",
    N=128,
    K=64,
    construction="pw", crc_poly="0x1D5",
)"""
"""CONFIG = dict(
    code_type="ps-pac",
    N=64,
    K=32,
    construction="dega",
    design_snr_db=2.0,
    alpha=7,
    pac_poly="1011011111",
)"""
# -----------------------------------------------------------------------

_FIVE_G = """0 1 2 4 8 16 32 3 5 64 9 6 17 10 18 128 12 33 65 20 256 34 24 36 7 129 66 512 11 40 68 130 19 13 48 14 72 257 21 132 35 258 26 513 80 37 25 22 136 260 264 38 514 96 67 41 144 28 69 42 516 49 74 272 160 520 288 528 192 544 70 44 131 81 50 73 15 320 133 52 23 134 384 76 137 82 56 27 97 39 259 84 138 145 261 29 43 98 515 88 140 30 146 71 262 265 161 576 45 100 640 51 148 46 75 266 273 517 104 162 53 193 152 77 164 768 268 274 518 54 83 57 521 112 135 78 289 194 85 276 522 58 168 139 99 86 60 280 89 290 529 524 196 141 101 147 176 142 530 321 31 200 90 545 292 322 532 263 149 102 105 304 296 163 92 47 267 385 546 324 208 386 150 153 165 106 55 328 536 577 548 113 154 79 269 108 578 224 166 519 552 195 270 641 523 275 580 291 59 169 560 114 277 156 87 197 116 170 61 531 525 642 281 278 526 177 293 388 91 584 769 198 172 120 201 336 62 282 143 103 178 294 93 644 202 592 323 392 297 770 107 180 151 209 284 648 94 204 298 400 608 352 325 533 155 210 305 547 300 109 184 534 537 115 167 225 326 306 772 157 656 329 110 117 212 171 776 330 226 549 538 387 308 216 416 271 279 158 337 550 672 118 332 579 540 389 173 121 553 199 784 179 228 338 312 704 390 174 554 581 393 283 122 448 353 561 203 63 340 394 527 582 556 181 295 285 232 124 205 182 643 562 286 585 299 354 211 401 185 396 344 586 645 593 535 240 206 95 327 564 800 402 356 307 301 417 213 568 832 588 186 646 404 227 896 594 418 302 649 771 360 539 111 331 214 309 188 449 217 408 609 596 551 650 229 159 420 310 541 773 610 657 333 119 600 339 218 368 652 230 391 313 450 542 334 233 555 774 175 123 658 612 341 777 220 314 424 395 673 583 355 287 183 234 125 557 660 616 342 316 241 778 563 345 452 397 403 207 674 558 785 432 357 187 236 664 624 587 780 705 126 242 565 398 346 456 358 405 303 569 244 595 189 566 676 361 706 589 215 786 647 348 419 406 464 680 801 362 590 409 570 788 597 572 219 311 708 598 601 651 421 792 802 611 602 410 231 688 653 248 369 190 364 654 659 335 480 315 221 370 613 422 425 451 614 543 235 412 343 372 775 317 222 426 453 237 559 833 804 712 834 661 808 779 617 604 433 720 816 836 347 897 243 662 454 318 675 618 898 781 376 428 665 736 567 840 625 238 359 457 399 787 591 678 434 677 349 245 458 666 620 363 127 191 782 407 436 626 571 465 681 246 707 350 599 668 790 460 249 682 573 411 803 789 709 365 440 628 689 374 423 466 793 250 371 481 574 413 603 366 468 655 900 805 615 684 710 429 794 252 373 605 848 690 713 632 482 806 427 904 414 223 663 692 835 619 472 455 796 809 714 721 837 716 864 810 606 912 722 696 377 435 817 319 621 812 484 430 838 667 488 239 378 459 622 627 437 380 818 461 496 669 679 724 841 629 351 467 438 737 251 462 442 441 469 247 683 842 738 899 670 783 849 820 728 928 791 367 901 630 685 844 633 711 253 691 824 902 686 740 850 375 444 470 483 415 485 905 795 473 634 744 852 960 865 693 797 906 715 807 474 636 694 254 717 575 913 798 811 379 697 431 607 489 866 723 486 908 718 813 476 856 839 725 698 914 752 868 819 814 439 929 490 623 671 739 916 463 843 381 497 930 821 726 961 872 492 631 729 700 443 741 845 920 382 822 851 730 498 880 742 445 471 635 932 687 903 825 500 846 745 826 732 446 962 936 475 853 867 637 907 487 695 746 828 753 854 857 504 799 255 964 909 719 477 915 638 748 944 869 491 699 754 858 478 968 383 910 815 976 870 917 727 493 873 701 931 756 860 499 731 823 922 874 918 502 933 743 760 881 494 702 921 501 876 847 992 447 733 827 934 882 937 963 747 505 855 924 734 829 965 938 884 506 749 945 966 755 859 940 830 911 871 639 888 479 946 750 969 508 861 757 970 919 875 862 758 948 977 923 972 761 877 952 495 703 935 978 883 762 503 925 878 735 993 885 939 994 980 926 764 941 967 886 831 947 507 889 984 751 942 996 971 890 509 949 973 1000 892 950 863 759 1008 510 979 953 763 974 954 879 981 982 927 995 765 956 887 985 997 986 943 891 998 766 511 988 1001 951 1002 893 975 894 1009 955 1004 1010 957 983 958 987 1012 999 1016 767 989 1003 990 1005 959 1011 1013 895 1006 1014 1017 1018 991 1020 1007 1015 1019 1021 1022 1023"""

def parse_index_set(text: str) -> list[int]:
    """Parse strings such as ``15,23,26-31,38,39``."""
    out: set[int] = set()
    for token in re.split(r"[\s,]+", text.strip()):
        if not token:
            continue
        if "-" in token:
            a, b = token.split("-", 1)
            lo, hi = int(a), int(b)
            if lo > hi:
                raise ValueError(f"bad range: {token}")
            out.update(range(lo, hi + 1))
        else:
            out.add(int(token))
    return sorted(out)


def parse_crc_polynomial(text: str) -> int:
    """Return q(x) as an integer bit mask; the MSB is the leading term."""
    s = text.strip().lower().replace(" ", "")
    if s.startswith(("0x", "0b", "0o")):
        q = int(s, 0)
    elif s.isdigit():
        q = int(s)
    else:
        q = 0
        for term in s.split("+"):
            if term == "1":
                e = 0
            elif term == "x":
                e = 1
            else:
                m = re.fullmatch(r"x\^?(\d+)", term)
                if not m:
                    raise ValueError(f"cannot parse CRC term '{term}'")
                e = int(m.group(1))
            q ^= 1 << e
    if q <= 0 or q.bit_length() < 2 or (q & 1) == 0:
        raise ValueError("CRC polynomial must have degree >= 1 and a nonzero constant term")
    return q


def polynomial_string(q: int) -> str:
    terms: list[str] = []
    for e in range(q.bit_length() - 1, -1, -1):
        if not ((q >> e) & 1):
            continue
        terms.append("1" if e == 0 else ("x" if e == 1 else f"x^{e}"))
    return "+".join(terms)


def crc_remainder(data_bits: Sequence[int], q: int) -> list[int]:
    """Systematic CRC remainder r(x)=d(x)x^t mod q(x).

    The bit order is the natural order of the data coordinates in I_crc,
    matching Sec. II-C of the manuscript where c=[d|r].
    """
    t = q.bit_length() - 1
    qbits = [(q >> e) & 1 for e in range(t, -1, -1)]
    work = list(data_bits) + [0] * t
    for p in range(len(data_bits)):
        if work[p]:
            for j, bit in enumerate(qbits):
                work[p + j] ^= bit
    return work[-t:]


def row_weight(i: int) -> int:
    # For G_N=G_2^{\otimes n}: w(g_i)=2^{wt(bin(i))}; manuscript Sec. II-B.
    return 1 << i.bit_count()


def core_set(i: int, info: set[int]) -> list[int]:
    """K_i in the Rowshan--Dau--Viterbo MWC formation formula."""
    return sorted(j for j in info if j > i and (j & ~i).bit_count() == 1)


def balancing_set(i: int, J: Sequence[int]) -> set[int]:
    """Construct M(J) from the minimum-weight formation theorem over F_2.

    This is the transparent leaf-level construction.  Algorithm 1 in the
    manuscript can update M incrementally during DFS; the candidate family and
    CRC test are the same.
    """
    groups: dict[int, list[int]] = {}
    for j in J:
        extra = j & ~i
        if extra == 0 or extra & (extra - 1):
            raise ValueError("J contains an index outside K_i")
        groups.setdefault(extra.bit_length() - 1, []).append(j)

    keys = sorted(groups)
    parity: dict[int, int] = {}

    def visit(pos: int, chosen: int, common: int, extra_or: int) -> None:
        if pos == len(keys):
            if chosen >= 2:
                m = extra_or | common
                parity[m] = parity.get(m, 0) ^ 1
            return
        visit(pos + 1, chosen, common, extra_or)
        b = keys[pos]
        for j in groups[b]:
            visit(pos + 1, chosen + 1, common & j, extra_or | (1 << b))

    visit(0, 0, i, 0)
    return {m for m, bit in parity.items() if bit}


@dataclass
class PolarCount:
    w_min: int
    multiplicity: int
    per_coset: dict[int, int]
    leaders: list[int]


def polar_mwc_count(info_indices: Iterable[int]) -> PolarCount:
    """A_wmin=sum_{i in B} 2^{|K_i|} (Rowshan--Dau--Viterbo)."""
    info = set(info_indices)
    if not info:
        raise ValueError("empty information set")
    w_min = min(row_weight(i) for i in info)
    leaders = sorted(i for i in info if row_weight(i) == w_min)
    per = {i: 1 << len(core_set(i, info)) for i in leaders}
    return PolarCount(w_min, sum(per.values()), per, leaders)


def polar_count_at_row_weight(info_indices: Iterable[int], weight: int) -> PolarCount:
    """Same coset formula restricted to leaders whose row weight is ``weight``."""
    info = set(info_indices)
    leaders = sorted(i for i in info if row_weight(i) == weight)
    per = {i: 1 << len(core_set(i, info)) for i in leaders}
    return PolarCount(weight, sum(per.values()), per, leaders)


def _dominates(i: int, j: int) -> bool:
    """Return True if the partial order makes j no worse than i."""
    if i == j or (i & ~j) == 0:
        return True
    if i.bit_count() != j.bit_count():
        return False
    width = max(i.bit_length(), j.bit_length()) + 1
    si = [b for b in range(width) if (i >> b) & 1]
    sj = [b for b in range(width) if (j >> b) & 1]
    return all(a <= b for a, b in zip(si, sj))


def pop_violations(N: int, info_indices: Iterable[int], limit: int = 8) -> list[tuple[int, int]]:
    """Examples (i,j) where selected i has a structurally better frozen j."""
    info = set(info_indices)
    frozen = set(range(N)) - info
    bad: list[tuple[int, int]] = []
    for i in sorted(info):
        for j in sorted(frozen):
            if _dominates(i, j):
                bad.append((i, j))
                if len(bad) >= limit:
                    return bad
    return bad


@dataclass
class CRCCount:
    reference_weight: int
    multiplicity_at_reference_weight: int
    per_coset: dict[int, int]
    leaders: list[int]
    candidates_tested: int
    candidates_rejected_by_frozen_rows: int


def crc_polar_mwc_count(
    N: int,
    info_indices: Iterable[int],
    q: int,
    crc_positions: Iterable[int] | None = None,
) -> CRCCount:
    """Count CRC-polar codewords at the underlying polar minimum row weight.

    Sec. III / Algorithm 1 of the accepted manuscript: enumerate the unique
    polar MWC combinations O={i} U J U M(J), discard combinations requiring a
    frozen row, and retain those whose values at R equal the systematic CRC
    remainder.  No message-by-message codeword enumeration is used.
    """
    info = set(info_indices)
    if not info:
        raise ValueError("empty CRC-polar information set")
    if min(info) < 0 or max(info) >= N:
        raise ValueError("information-set index outside [0,N-1]")

    t = q.bit_length() - 1
    R = set(sorted(info)[-t:] if crc_positions is None else crc_positions)
    if len(R) != t or not R.issubset(info):
        raise ValueError(f"CRC position set must contain exactly t={t} positions inside I_crc")

    data = sorted(info - R)
    if data and R and max(data) >= min(R):
        raise ValueError("requires the manuscript placement max(data)<min(R)")
    R_order = sorted(R)

    w_min = min(row_weight(i) for i in info)
    leaders = sorted(i for i in info if row_weight(i) == w_min)
    per: dict[int, int] = {}
    tested = 0
    rejected_frozen = 0

    for i in leaders:
        # Sec. III: a CRC-constrained coordinate cannot be a free coset leader.
        if i in R:
            per[i] = 0
            continue

        Ki = core_set(i, info)
        count_i = 0
        for mask in range(1 << len(Ki)):
            J = [Ki[b] for b in range(len(Ki)) if (mask >> b) & 1]
            M = balancing_set(i, J)
            O = {i, *J, *M}
            tested += 1

            if not O.issubset(info):
                rejected_frozen += 1
                continue

            rem = crc_remainder([int(d in O) for d in data], q)
            required = [int(r in O) for r in R_order]
            if rem == required:
                count_i += 1
        per[i] = count_i

    return CRCCount(w_min, sum(per.values()), per, leaders, tested, rejected_frozen)


def _phi(x: float) -> float:
    if x <= 0.0:
        return 1.0
    if x < 10.0:
        return math.exp(-0.4527 * (x**0.86) + 0.0218)
    v = math.sqrt(math.pi / x) * (1.0 - 10.0 / (7.0 * x)) * math.exp(-x / 4.0)
    return max(0.0, min(1.0, v))


def _phi_inv(y: float) -> float:
    y = min(max(y, 1e-300), 1.0 - 1e-15)
    lo, hi = 0.0, 2048.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _phi(mid) > y:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def dega_reliabilities(N: int, design_snr_db: float, design_rate: float) -> list[float]:
    """Gaussian-approximation mean LLRs for a BI-AWGN polar construction."""
    if N <= 0 or N & (N - 1):
        raise ValueError("N must be a power of two")
    if not (0.0 < design_rate <= 1.0):
        raise ValueError("design rate must be in (0,1]")

    mean0 = 4.0 * design_rate * (10.0 ** (design_snr_db / 10.0))
    means = [mean0]
    for _ in range(int(math.log2(N))):
        nxt: list[float] = []
        for m in means:
            p = _phi(m)
            nxt.extend((_phi_inv(1.0 - (1.0 - p) ** 2), 2.0 * m))
        means = nxt
    return means


def reliability_order_dega(N: int, design_snr_db: float, design_rate: float) -> list[int]:
    means = dega_reliabilities(N, design_snr_db, design_rate)
    return sorted(range(N), key=lambda i: (means[i], i))


def construct_dega(N: int, size: int, design_snr_db: float, design_rate: float) -> list[int]:
    return sorted(reliability_order_dega(N, design_snr_db, design_rate)[-size:])


def pw_reliabilities(N: int, beta: float = 2.0 ** 0.25) -> list[float]:
    """Polarization-weight scores; larger score means more reliable.

    PW supplies a separate reliability-based information-set construction.
    It is not one of the three main constructors used by this repository.
    """
    if N <= 0 or N & (N - 1):
        raise ValueError("N must be a power of two")
    n = int(math.log2(N))
    return [sum(((i >> j) & 1) * (beta**j) for j in range(n)) for i in range(N)]


def reliability_order_pw(N: int, beta: float = 2.0 ** 0.25) -> list[int]:
    scores = pw_reliabilities(N, beta)
    return sorted(range(N), key=lambda i: (scores[i], i))


def construct_pw(N: int, size: int, beta: float = 2.0 ** 0.25) -> list[int]:
    return sorted(reliability_order_pw(N, beta)[-size:])


def construct_crossed_rm(N: int, size: int, reliability_order: Sequence[int]) -> list[int]:
    """Crossed RM-polar profile of Rowshan--Burg--Viterbo (2021).

    Let r be the largest integer for which sum_{j=0}^r C(n,j) <= K.  All
    generator rows with weight at least 2^(n-r) are selected.  If fewer than K
    rows have been selected, the remaining positions are the most reliable
    rows from the next RM weight class.  ``reliability_order`` is least to most
    reliable; DEGA is used by the CLI to reproduce the 2021 construction.
    """
    if N <= 0 or N & (N - 1):
        raise ValueError("N must be a power of two")
    if not (0 < size <= N):
        raise ValueError("profile size must satisfy 0<size<=N")
    if len(reliability_order) != N or set(reliability_order) != set(range(N)):
        raise ValueError("reliability_order must be a permutation of 0..N-1")

    n = int(math.log2(N))
    cumulative = 0
    r = -1
    for candidate in range(n + 1):
        trial = cumulative + math.comb(n, candidate)
        if trial > size:
            break
        cumulative = trial
        r = candidate

    min_popcount = n - r
    chosen = {i for i in range(N) if i.bit_count() >= min_popcount}
    need = size - len(chosen)
    if need:
        boundary_popcount = min_popcount - 1
        rank = {i: pos for pos, i in enumerate(reliability_order)}
        boundary = [i for i in range(N) if i.bit_count() == boundary_popcount]
        boundary.sort(key=lambda i: rank[i], reverse=True)
        chosen.update(boundary[:need])

    if len(chosen) != size:
        raise RuntimeError("crossed RM-polar construction produced the wrong dimension")
    return sorted(chosen)

def reliability_order_5g(N: int) -> list[int]:
    if N > 1024:
        raise ValueError("5G universal sequence requires N<=1024")
    order = [int(s) for s in _FIVE_G.split() if int(s)<N]
    if len(order)!=N or set(order)!=set(range(N)):
        raise ValueError("invalid embedded 5G reliability sequence")
    return order


def profile_from_construction(name, N, size, *, snr=2.0, rate=0.5, beta=2**0.25):
    if not (0 < size <= N):
        raise ValueError("profile size must be in [1,N]")
    if name=="5g":
        order = reliability_order_5g(N)
        return sorted(order[-size:]), order
    if name=="dega":
        order = reliability_order_dega(N,snr,rate)
        return sorted(order[-size:]), order
    if name=="crossed-rm":
        order = reliability_order_dega(N,snr,rate)
        return construct_crossed_rm(N,size,order), order
    if name=="pw":
        order = reliability_order_pw(N,beta)
        return sorted(order[-size:]), order
    raise ValueError("unknown construction")

def parse_convolution_polynomial(text: str) -> list[int]:
    """Parse p=[p0,p1,...,ps], with p0 multiplying the current bit.

    Accepted forms include ``1011011``, ``0b1011011`` and
    ``[1,0,1,1,0,1,1]``.  The left-most bit is p0, matching

        u_i = sum_{ell=0}^s p_ell v_{i-ell}.
    """
    s = text.strip().lower()
    if s.startswith("0b"):
        s = s[2:]
    if re.fullmatch(r"[01]+", s):
        bits = [int(c) for c in s]
    else:
        bits = [int(x) for x in re.findall(r"[01]", s)]
    if not bits:
        raise ValueError("cannot parse convolution polynomial")
    if bits[0] != 1:
        raise ValueError("PAC precoder requires p0=1")
    if bits[-1] != 1:
        raise ValueError("the last coefficient should be nonzero; trim trailing zeros")
    return bits


def polynomial_bits(p: Sequence[int]) -> str:
    return "".join(str(x) for x in p)


def shifted_information_set(
    base_info: Iterable[int],
    reliability_order: Sequence[int],
    alpha: int,
) -> tuple[list[int], list[int], list[int]]:
    """Apply the PS-PAC profile shift from the TCOM manuscript.

    ``reliability_order`` is least -> most reliable.  We freeze the ``alpha``
    most reliable coordinates currently in the information set and enable the
    ``alpha`` most reliable coordinates of the old frozen set.  For a standard
    top-K reliability profile, this is exactly the shifted window in Sec. IV.
    The same operation is also useful with a crossed RM-polar base profile.

    Returns ``(I_shifted, removed_from_I, added_to_I)``.
    """
    order = list(reliability_order)
    N = len(order)
    base = set(base_info)
    if set(order) != set(range(N)):
        raise ValueError("reliability_order must be a permutation of 0..N-1")
    if not (0 <= alpha <= min(len(base), N - len(base))):
        raise ValueError("alpha is outside the valid profile-shift range")
    if alpha == 0:
        return sorted(base), [], []

    rank = {idx: pos for pos, idx in enumerate(order)}
    removed = sorted(base, key=lambda x: rank[x], reverse=True)[:alpha]
    frozen = set(range(N)) - base
    added = sorted(frozen, key=lambda x: rank[x], reverse=True)[:alpha]
    shifted = (base - set(removed)) | set(added)
    return sorted(shifted), sorted(removed), sorted(added)


def _add_to_m(i: int, j: int, J: list[int], M: set[int], O: bytearray) -> None:
    """Incremental update of the balancing set M after adding core row j.

    This is the ``addToM`` step used in the fast PAC enumerator.  Only pairs
    involving the newly added row can create a new balancing index.  Repeated
    generation of the same index cancels over F_2.
    """
    for source in (list(M), list(J)):
        for y in source:
            # j and y must introduce distinct coordinates outside supp(i).
            if ((j & y) & ~i) != 0:
                continue
            m = ((j | y) & ~i) | (i & j & y)
            if O[m]:
                O[m] = 0
                M.discard(m)
            else:
                O[m] = 1
                M.add(m)


def _convolution_memory(v: Sequence[int], z: int, p: Sequence[int]) -> int:
    """Contribution to u_z from v_{z-1},...,v_{z-s}; p0*v_z excluded."""
    bit = 0
    upto = min(len(p) - 1, z)
    for ell in range(1, upto + 1):
        if p[ell]:
            bit ^= v[z - ell]
    return bit


@dataclass
class PACCosetCount:
    leader: int
    row_weight: int
    multiplicity: int
    core_size: int
    truncated_core_size: int
    capable_frozen_count: int
    subsets_tested: int
    subsets_rejected: int


def pac_coset_mwc_count(
    N: int,
    info: set[int],
    p: Sequence[int],
    i: int,
) -> PACCosetCount:
    """Count row-weight-w(g_i) words in PAC coset C_i structurally.

    Following the capable-coset analysis of Rowshan--Yuan, let

        F* = {f in I^c, f>i : |supp(f)\\supp(i)|>1}.

    If F* is empty, the precoder cannot reduce the polar count in this coset.
    Otherwise, only core coordinates before f_max=max(F*) are explicitly
    enumerated.  This replaces a 2^K codebook search by a much smaller
    2^|K_i^f| search local to the coset.
    """
    Ki = core_set(i, info)
    frozen = set(range(N)) - info
    Fstar = sorted(f for f in frozen if f > i and (f & ~i).bit_count() > 1)

    if not Fstar:
        return PACCosetCount(
            i,
            row_weight(i),
            1 << len(Ki),
            len(Ki),
            0,
            0,
            0,
            0,
        )

    fmax = Fstar[-1]
    Fstar_set = set(Fstar)
    Kif = [j for j in Ki if j < fmax]
    total_subsets = 1 << len(Kif)
    valid_subsets = 0
    rejected = 0

    for mask in range(total_subsets):
        chosen_core = {Kif[b] for b in range(len(Kif)) if (mask >> b) & 1}
        J: list[int] = []
        M: set[int] = set()
        O = bytearray(N)
        O[i] = 1

        # v is the rate-profile vector before convolution.  Positions before
        # the coset leader i are zero by definition of the coset.
        v = bytearray(N)
        v[i] = 1
        bad = False

        for z in range(i + 1, fmax + 1):
            memory = _convolution_memory(v, z, p)

            if z in info:
                if z in chosen_core:
                    _add_to_m(i, z, J, M, O)
                    J.append(z)
                    desired_u = 1
                else:
                    desired_u = O[z]

                # p0=1, hence v_z = u_z XOR memory.
                v[z] = memory ^ desired_u
            else:
                # v_z is frozen to zero, but convolution may make u_z=1.
                u_z = memory
                if z in Fstar_set:
                    if u_z != O[z]:
                        rejected += 1
                        bad = True
                        break
                elif u_z:
                    # A frozen core-type coordinate forced to one becomes part
                    # of the row combination and must update M incrementally.
                    _add_to_m(i, z, J, M, O)
                    J.append(z)

        if not bad:
            valid_subsets += 1

    multiplier = 1 << (len(Ki) - len(Kif))
    return PACCosetCount(
        leader=i,
        row_weight=row_weight(i),
        multiplicity=valid_subsets * multiplier,
        core_size=len(Ki),
        truncated_core_size=len(Kif),
        capable_frozen_count=len(Fstar),
        subsets_tested=total_subsets,
        subsets_rejected=rejected,
    )


@dataclass
class PACWeightCount:
    weight: int
    multiplicity: int
    per_coset: dict[int, int]
    leaders: list[int]
    subsets_tested: int
    subsets_rejected: int


def pac_count_at_weight(
    N: int,
    info_indices: Iterable[int],
    p: Sequence[int],
    weight: int,
) -> PACWeightCount:
    """Count the PAC minimum-row-weight class for leaders of given weight.

    For ``weight=min_i w(g_i)``, this is the exact A_wmin enumerated by the
    fast PAC method.  Supplying a larger row-weight is useful for reproducing
    manuscript tables when the lower row-weight class is known to be empty;
    it should not be interpreted as a general higher-weight spectrum routine.
    """
    info = set(info_indices)
    leaders = sorted(i for i in info if row_weight(i) == weight)
    per: dict[int, int] = {}
    tested = 0
    rejected = 0
    for i in leaders:
        r = pac_coset_mwc_count(N, info, p, i)
        per[i] = r.multiplicity
        tested += r.subsets_tested
        rejected += r.subsets_rejected
    return PACWeightCount(weight, sum(per.values()), per, leaders, tested, rejected)

# ============================================================================
# Exact low-weight spectrum supplement: parity-check syndrome meet-in-the-middle
# ============================================================================
# This searches coordinate supports (not data messages / the full 2^K codebook).
# Once structural counting eliminates the row-weight lower bound, matching
# syndromes of disjoint coordinate halves counts every codeword of weight w
# exactly once. It can be exponential in w; explicit limits prevent runaway.


def polar_row_bits(i: int) -> int:
    """Bitset of row i in F^{tensor log2(N)} (natural order)."""
    word = 0
    j = i
    while True:
        word |= 1 << j
        if j == 0:
            break
        j = (j - 1) & i
    return word


def binary_generator(code_type: str, N: int, info: list[int],
                     poly: list[int] | None = None,
                     crc_q: int | None = None,
                     crc_positions: list[int] | None = None) -> list[int]:
    """Actual code's K generator rows as Python integer bitsets."""
    row_g = [polar_row_bits(i) for i in range(N)]
    if code_type == 'polar':
        return [row_g[i] for i in info]
    if code_type in ('pac', 'ps-pac'):
        if not poly:
            raise ValueError('missing convolutional polynomial')
        return [
            _xor_pac_generator(row_g, i, N, poly)
            for i in info
        ]
    if code_type == 'crc-polar':
        if crc_q is None or crc_positions is None:
            raise ValueError('CRC polynomial/positions are required')
        R = sorted(crc_positions)
        data = sorted(set(info) - set(R))
        if len(R) != crc_q.bit_length()-1 or (data and max(data)>=min(R)):
            raise ValueError('CRC bit placement must follow systematic data positions')
        out = []
        for d_pos, di in enumerate(data):
            basis = [0]*len(data)
            basis[d_pos] = 1
            rem = crc_remainder(basis, crc_q)
            row = row_g[di]
            for j,bit in enumerate(rem):
                if bit:
                    row ^= row_g[R[j]]
            out.append(row)
        return out
    raise ValueError(f'unsupported code type: {code_type}')


def _xor_pac_generator(rows: list[int], i: int, N: int, p: list[int]) -> int:
    out = 0
    for ell, bit in enumerate(p):
        if i + ell >= N:
            break
        if bit:
            out ^= rows[i + ell]
    return out


def parity_check_columns(gen: list[int], N: int) -> list[int]:
    """RREF-derived H columns; x belongs to C iff XOR_{j:x_j=1} H_j=0."""
    K = len(gen)
    matrix = list(gen)
    pivots = []
    rank = 0
    for j in range(N):
        row = next((r for r in range(rank,K) if ((matrix[r]>>j)&1)),None)
        if row is None:
            continue
        matrix[rank],matrix[row] = matrix[row],matrix[rank]
        pivot = matrix[rank]
        for r in range(K):
            if r != rank and ((matrix[r] >> j)&1):
                matrix[r] ^= pivot
        pivots.append(j)
        rank += 1
        if rank == K:
            break
    if rank != K:
        raise ValueError(f'code generator has rank {rank} not {K}; check parameters')
    pivot_set = set(pivots)
    free = [j for j in range(N) if j not in pivot_set]
    cols = [0]*N
    for f_idx, f in enumerate(free):
        v = 1 << f_idx
        cols[f] = v
        for r, pi in enumerate(pivots):
            if (matrix[r] >> f)&1:
                cols[pi] |= v
    return cols


def _iter_syndromes(columns: Sequence[int], take: int):
    """Generate support syndromes for fixed-weight selections in one half.

    itertools.combinations executes the combinatorial iterator in CPython's
    C implementation. Each yielded XOR is computed with Python integers,
    so parity-check matrices with more than 64 rows are supported.
    """
    if take == 0:
        yield 0
        return
    for combo in itertools.combinations(columns, take):
        syndrome = 0
        for col in combo:
            syndrome ^= col
        yield syndrome


def _weight_syndrome_count(columns: Sequence[int], weight: int,
                           max_combos: int, max_hash: int) -> int:
    """Count exact weight-w codewords by disjoint parity-check syndrome matches.

    The code coordinate set is partitioned into two fixed halves. For each
    split l+r=w, a dictionary counts syndromes on the smaller half and is
    queried using combinations on the other half. Each support is counted
    precisely once, even when multiple supports share a syndrome.
    """
    N = len(columns)
    left_len = N // 2
    right_len = N - left_len
    left = columns[:left_len]
    right = columns[left_len:]
    plan = []
    for l in range(max(0, weight - right_len), min(weight, left_len) + 1):
        r = weight - l
        nl = math.comb(left_len, l)
        nr = math.comb(right_len, r)
        if max(nl, nr) > max_combos:
            raise ResourceLimit(
                f'weight {weight}: split {l}+{r} needs '
                f'{max(nl, nr):,} support combinations (limit {max_combos:,})'
            )
        if min(nl, nr) > max_hash:
            raise ResourceLimit(
                f'weight {weight}: split {l}+{r} needs '
                f'{min(nl, nr):,} hash entries (limit {max_hash:,})'
            )
        plan.append((l, r, nl, nr))
    answer = 0
    for l, r, nl, nr in plan:
        if nl <= nr:
            small, small_weight = left, l
            large, large_weight = right, r
        else:
            small, small_weight = right, r
            large, large_weight = left, l
        lookup: dict[int, int] = {}
        for syndrome in _iter_syndromes(small, small_weight):
            lookup[syndrome] = lookup.get(syndrome, 0) + 1
        for syndrome in _iter_syndromes(large, large_weight):
            answer += lookup.get(syndrome, 0)
    return answer


class ResourceLimit(RuntimeError):
    """A rigorous minimum-distance certificate exceeds configured resources."""


def exact_low_weight_search(gen: list[int], N: int, min_weight: int,
                            upper_weight: int, *, max_weight: int = 10,
                            max_combos: int = 12_000_000,
                            max_hash: int = 1_000_000):
    """Certify A_dmin by syndrome counting, without enumerating 2^K messages.

    Counts weights sequentially; does not skip unexplored weights. A code whose
    generator rows have even parity has an all-even weight spectrum, so odd
    weights can be skipped with a proof. Results contain all verified counts.
    """
    truncated = upper_weight > max_weight
    search_stop = min(upper_weight,max_weight)
    cols = parity_check_columns(gen, N)
    all_even = all(x.bit_count()%2==0 for x in gen)
    checked={}
    for w in range(min_weight,search_stop+1):
        if all_even and w%2:
            checked[w]=0
            continue
        try:
            a = _weight_syndrome_count(cols,w,max_combos,max_hash)
        except ResourceLimit as e:
            return dict(d_min=None,A_dmin=None, checked=checked,reason=str(e))
        checked[w]=a
        if a:
            return dict(d_min=w,A_dmin=a,checked=checked,reason=None)
    return dict(d_min=None,A_dmin=None,checked=checked,
                reason=(f'no codeword through weight {max_weight}; larger weights not checked'
                        if truncated else 'no word within certified upper bound; check implementation'))


def _representable_row_class(w: int) -> bool:
    return w > 0 and w&(w-1)==0


def enumerate_code(*, code_type: str, N: int, K: int, construction: str='dega',
                   info_set: str | Sequence[int] | None=None,
                   reliability_order: str | Sequence[int] | None=None,
                   alpha: int=0, crc_poly: str='x^5+x^3+1',
                   pac_poly: str='1011011', crc_positions: str | Sequence[int] | None=None,
                   design_snr_db:float=2.0, design_rate:float|None=None,
                   pw_beta:float=2**0.25,
                   max_exact_weight:int=10, max_syndrome_combos:int=12_000_000,
                   max_hash_entries:int=1_000_000,
                   max_structural_subsets:int=16_777_216,
                   per_coset: bool=False) -> dict:
    """Public API. Only one selected family is enumerated per call.

    Returns A_dmin as integer when *certified*, otherwise None with reason.
    Input K is data dimension for all four families (excludes CRC parity bits).
    Explicit `info_set` for PS-PAC means the shifted information set if no
    reliability order is supplied; otherwise it is the unshifted/base set.
    """
    start=time.perf_counter()
    typ=code_type.strip().lower().replace('_','-')
    if typ not in ('polar','pac','ps-pac','crc-polar'):
        raise ValueError('code_type must be polar, pac, ps-pac or crc-polar')
    if N < 2 or (N & (N-1)):
        raise ValueError('N must be a power of two')
    if not (0<K<N):
        raise ValueError('K must satisfy 0<K<N')
    q=parse_crc_polynomial(crc_poly) if typ=='crc-polar' else None
    t=q.bit_length()-1 if q else 0
    size=K+t
    if size > N:
        raise ValueError('K+CRC degree cannot exceed N')
    if design_rate is None:
        design_rate=K/N
    if construction=='explicit':
        if info_set is None:
            raise ValueError('explicit construction requires info_set')
        I=parse_index_set(info_set) if isinstance(info_set,str) else sorted(set(info_set))
        order=None
        if reliability_order is not None:
            order=(parse_index_set(reliability_order) if isinstance(reliability_order,str)
                   else list(reliability_order))
            # A reliability sequence is an ordered permutation; do not use
            # parse_index_set if it was given as a string (which sorts it).
            if isinstance(reliability_order,str):
                order=[int(x) for x in re.split(r'[\s,]+', reliability_order.strip()) if x]
    else:
        I, order=profile_from_construction(construction,N,size,snr=design_snr_db,
                                           rate=design_rate,beta=pw_beta)
        if info_set is not None:
            raise ValueError('info_set is only used with construction=explicit')
    removed=[]; added=[]
    if typ=='ps-pac':
        if construction=='explicit' and order is None:
            # Already-shifted profile supplied verbatim; alpha is informational.
            pass
        else:
            I, removed, added=shifted_information_set(I,order,alpha)
    if len(I)!=size or min(I)<0 or max(I)>=N:
        raise ValueError(f'information set must contain exactly {size} unique indices in [0,N)')
    if typ=='crc-polar':
        R=(parse_index_set(crc_positions) if isinstance(crc_positions,str) else
           sorted(set(crc_positions)) if crc_positions is not None else sorted(I)[-t:])
        if len(R)!=t or not set(R).issubset(I):
            raise ValueError('CRC positions must be t distinct selected indices')
        if max(set(I)-set(R))>=min(R):
            raise ValueError('CRC positions must follow data indices in natural order')
    else:
        R=None
    p=parse_convolution_polynomial(pac_poly) if typ in ('pac','ps-pac') else None
    pop_ok = not pop_violations(N,I)
    setup=time.perf_counter()-start
    t0=time.perf_counter()
    lower=min(row_weight(i) for i in I)
    structural=None; weights=[]; upper_candidate=None; cosine_counts={}
    if typ=='polar':
        vio=pop_violations(N,I)
        if vio:
            # Polar formula A_dmin for an arbitrary profile may not be exact.
            # Exact syndrome fallback for modest instances.
            structural_status='partial order violated; structural multiplicity unavailable'
        else:
            f=polar_mwc_count(I)
            return dict(code_type=typ,N=N,K=K,construction=construction,
                        info_set=I, d_min=f.w_min,A_dmin=f.multiplicity,
                        certified=True,method='Rowshan-Dau-Viterbo polar MWC formula',
                        lower_bound=f.w_min, per_coset=f.per_coset if per_coset else None,
                        weight_counts={f.w_min:f.multiplicity},
                        seconds=dict(profile=setup,enumeration=time.perf_counter()-t0))
    elif typ=='crc-polar':
        violated=pop_violations(N,I)
        if violated:
            structural_status='partial order violated; CRC structural enumeration unavailable'
        else:
            max_core=max([len(core_set(i,set(I))) for i in I if row_weight(i)==lower and i not in set(R)] or [0])
            if (1<<max_core)>max_structural_subsets:
                structural_status=f'CRC structural core would need 2^{max_core} subsets; cap reached'
            else:
                f=crc_polar_mwc_count(N,I,q,R)
                weights.append((lower,f.multiplicity_at_reference_weight))
                cosine_counts[lower]=f.per_coset
                if f.multiplicity_at_reference_weight:
                    return dict(code_type=typ,N=N,K=K,construction=construction,
                                info_set=I,crc_positions=R,crc_polynomial=polynomial_string(q),
                                d_min=lower,A_dmin=f.multiplicity_at_reference_weight,
                                certified=True,method='CRC-surviving polar MWC structural enumeration',
                                lower_bound=lower, per_coset=f.per_coset if per_coset else None,
                                weight_counts={lower:f.multiplicity_at_reference_weight},
                                seconds=dict(profile=setup,enumeration=time.perf_counter()-t0))
                structural_status='minimum row-weight class was eliminated by CRC'
    elif typ in ('pac','ps-pac'):
        structural_status=''
        for w in sorted(set(row_weight(i) for i in I)):
            sizes=[len([j for j in core_set(i,set(I)) if j < max([f for f in range(N) if f not in I and f>i and (f & ~i).bit_count()>1] or [0])])
                   for i in I if row_weight(i)==w]
            if sizes and ((1<<max(sizes)) > max_structural_subsets):
                structural_status=f'PAC structural core exceeds {max_structural_subsets:,} subset cap'
                break
            f=pac_count_at_weight(N,I,p,w)
            weights.append((w,f.multiplicity))
            cosine_counts[w]=f.per_coset
            if f.multiplicity:
                upper_candidate=w
                if w==lower and pop_ok:
                    return dict(code_type=typ,N=N,K=K,construction=construction,
                                info_set=I,alpha=alpha if typ=='ps-pac' else None,
                                removed=removed,added=added,pac_polynomial=polynomial_bits(p),
                                d_min=w,A_dmin=f.multiplicity,certified=True,
                                method='capable-coset structural MWC enumeration',lower_bound=lower,
                                weight_counts={w:f.multiplicity},
                                per_coset=f.per_coset if per_coset else None,
                                seconds=dict(profile=setup,enumeration=time.perf_counter()-t0))
                structural_status=('lower-weight structural class eliminated; exact weight search required'
                    if w>lower else 'non-POP profile: independently checking the minimum-weight multiplicity')
                break
    # When the structural reference weight is removed, derive an ordinary
    # binary generator and *exactly* count supports, not data messages.
    gen=binary_generator(typ,N,I,p,q,R)
    upper=min(x.bit_count() for x in gen)
    if upper_candidate is not None:
        upper=min(upper,upper_candidate)
    if typ=='polar' and not weights:
        structural_status='non-POP polar profile; exact low-weight search required'
    exact=exact_low_weight_search(gen,N,lower,upper,max_weight=max_exact_weight,
                                  max_combos=max_syndrome_combos,max_hash=max_hash_entries)
    certified=exact['d_min'] is not None
    return dict(code_type=typ,N=N,K=K,construction=construction,info_set=I,
                alpha=alpha if typ=='ps-pac' else None,
                removed=removed,added=added,
                pac_polynomial=polynomial_bits(p) if p else None,
                crc_polynomial=polynomial_string(q) if q else None,
                crc_positions=R,d_min=exact['d_min'],A_dmin=exact['A_dmin'],
                certified=certified,method='exact syndrome support enumeration' if certified else 'unresolved',
                lower_bound=lower,upper_bound=upper,
                weight_counts=exact['checked'],structural_classes=dict(weights),
                status=(structural_status if certified else (structural_status+'; '+str(exact['reason']))),
                per_coset=cosine_counts if per_coset else None,
                seconds=dict(profile=setup,enumeration=time.perf_counter()-t0))


def _arg_parser():
    ap=argparse.ArgumentParser(description='Polar / PAC / PS-PAC / CRC-polar minimum distance and multiplicity')
    ap.add_argument('--code-type',choices=['polar','pac','ps-pac','crc-polar'],required=True)
    ap.add_argument('--N',type=int,required=True)
    ap.add_argument('--K',type=int,required=True,help='information bits (excludes CRC parity bits)')
    ap.add_argument('--construction',choices=['dega','5g','crossed-rm','pw','explicit'],default='dega')
    ap.add_argument('--info-set',help='explicit sorted indices/ranges, e.g. 15,23,26-31')
    ap.add_argument('--reliability-order',help='least-to-most reliable permutation for explicit PS-PAC shift')
    ap.add_argument('--alpha',type=int,default=7,help='PS-PAC shift (ignored for other code types)')
    ap.add_argument('--pac-poly',default='1011011111',help='PAC/PS-PAC convolution coefficients p0,..,ps')
    ap.add_argument('--crc-poly',default='x^5+x^3+1')
    ap.add_argument('--crc-positions',help='CRC indices R, default last t indices in I_crc')
    ap.add_argument('--design-snr-db',type=float,default=2.0)
    ap.add_argument('--design-rate',type=float,help='default K/N')
    ap.add_argument('--pw-beta',type=float,default=2**0.25)
    ap.add_argument('--max-exact-weight',type=int,default=10)
    ap.add_argument('--max-syndrome-combos',type=int,default=12_000_000)
    ap.add_argument('--max-hash-entries',type=int,default=1_000_000)
    ap.add_argument('--max-structural-subsets',type=int,default=16_777_216)
    ap.add_argument('--per-coset',action='store_true')
    return ap


def _print_result(r):
    print(f"Code: {r['code_type']}  N={r['N']}  K={r['K']}  profile={r['construction']}")
    print(f"Information set ({len(r['info_set'])}): {r['info_set']}")
    if r.get('alpha') is not None:
        print(f"Shift alpha={r['alpha']}  removed={r['removed']}  added={r['added']}")
    if r.get('pac_polynomial'):
        print('Convolution polynomial p:',r['pac_polynomial'])
    if r.get('crc_polynomial'):
        print('CRC polynomial q:',r['crc_polynomial'])
    print('Polar-coset lower bound:',r['lower_bound'])
    if r.get('structural_classes'):
        print('Structural row-weight contributions:',r['structural_classes'])
    print('Verified weight counts:',r.get('weight_counts',{}))
    if r['certified']:
        print(f"d_min = {r['d_min']}")
        print(f"A_dmin = {r['A_dmin']}")
        print('CERTIFIED using:',r['method'])
    else:
        print('d_min = NOT CERTIFIED; A_dmin = NOT CERTIFIED')
        print('status:',r.get('status','resource or method limitation'))
        if r.get('upper_bound'):
            print(f"known distance bounds: {r['lower_bound']} <= d_min <= {r['upper_bound']}")
    if r.get('per_coset') is not None:
        print('Per-coset:',r['per_coset'])
    print(f"Elapsed enumeration: {r['seconds']['enumeration']:.4f} s")


def main(argv=None):
    """Run CONFIG with no arguments; any arguments select the CLI instead."""
    if argv is None:
        argv = sys.argv[1:]
    if not argv and CONFIG is not None:
        result = enumerate_code(**CONFIG)
        _print_result(result)
        return result
    args = vars(_arg_parser().parse_args(argv))
    result = enumerate_code(**args)
    _print_result(result)
    return result

if __name__=='__main__':
    main()
