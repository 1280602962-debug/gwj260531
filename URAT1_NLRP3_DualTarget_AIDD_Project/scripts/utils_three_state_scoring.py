#!/usr/bin/env python3
"""Boltzmann pi weights and S_pi for URAT1 three-state benchmark (Protocol C)."""

from __future__ import annotations

import math
from typing import Mapping

# kcal/mol at 298 K
RT_KCAL_MOL = 0.593

STATE_INWARD = "inward"
STATE_OCCLUDED = "occluded"
STATE_OUTWARD = "outward"
DEFAULT_STATES = (STATE_INWARD, STATE_OCCLUDED, STATE_OUTWARD)


def boltzmann_pi(
    delta_g: Mapping[str, float],
    temperature_rt: float = RT_KCAL_MOL,
) -> dict[str, float]:
    """Compute Boltzmann state populations from rescored free energies (kcal/mol)."""
    states = [s for s in DEFAULT_STATES if s in delta_g]
    if not states:
        raise ValueError("delta_g must include at least one of inward/occluded/outward")
    dg = [float(delta_g[s]) for s in states]
    dg_min = min(dg)
    weights = [math.exp(-(g - dg_min) / temperature_rt) for g in dg]
    total = sum(weights)
    if total <= 0.0:
        uniform = 1.0 / len(states)
        return {s: uniform for s in states}
    return {s: w / total for s, w in zip(states, weights)}


def s_pi(pi: Mapping[str, float]) -> float:
    """Directional three-state score: favor inward + occluded over outward."""
    return float(pi.get(STATE_INWARD, 0.0) + pi.get(STATE_OCCLUDED, 0.0) - pi.get(STATE_OUTWARD, 0.0))


def apply_outward_clash_penalty(
    dg_outward: float,
    clash: bool,
    penalty_kcal: float = 5.0,
) -> float:
    """Add fixed penalty when copied outward pose has severe vdW clash."""
    return dg_outward + (penalty_kcal if clash else 0.0)


def score_molecule(
    dg_inward: float,
    dg_occluded: float,
    dg_outward: float,
    outward_clash: bool = False,
    outward_penalty_kcal: float = 5.0,
    temperature_rt: float = RT_KCAL_MOL,
) -> dict[str, float | bool]:
    """Full Protocol C post-processing for one compound."""
    dg_out = apply_outward_clash_penalty(dg_outward, outward_clash, outward_penalty_kcal)
    delta_g = {
        STATE_INWARD: dg_inward,
        STATE_OCCLUDED: dg_occluded,
        STATE_OUTWARD: dg_out,
    }
    pi = boltzmann_pi(delta_g, temperature_rt=temperature_rt)
    return {
        "dg_inward": dg_inward,
        "dg_occluded": dg_occluded,
        "dg_outward": dg_out,
        "pi_inward": pi[STATE_INWARD],
        "pi_occluded": pi[STATE_OCCLUDED],
        "pi_outward": pi[STATE_OUTWARD],
        "s_pi": s_pi(pi),
        "outward_clash": outward_clash,
        "penalty_applied": outward_penalty_kcal if outward_clash else 0.0,
    }


if __name__ == "__main__":
    # Example: inward-favoring inhibitor-like profile
    example = score_molecule(
        dg_inward=-45.0,
        dg_occluded=-42.0,
        dg_outward=-38.0,
        outward_clash=True,
    )
    for key, value in example.items():
        print(f"{key}: {value}")
