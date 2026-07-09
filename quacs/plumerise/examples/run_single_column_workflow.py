"""Run one synthetic profile through the wildfire profile fraction workflow."""

import json

import numpy as np

from quacs.plumerise import compute_wildfire_profile_fraction_driver


def build_synthetic_profile_inputs():
    z = np.arange(0.0, 16000.0, 1000.0)
    t = np.where(z <= 10000.0, 300.0 - 6.5 * (z / 1000.0), 235.0)
    p = 1013.25 * np.exp(-z / 8000.0)
    return {
        "z": z,
        "p": p,
        "t": t,
        "u": np.ones_like(z),
        "v": np.zeros_like(z),
        "qv": np.full_like(z, 0.01),
    }


def main():
    layer_fraction = compute_wildfire_profile_fraction_driver(
        **build_synthetic_profile_inputs(),
        vegetation_class="forest",
        fire_size_mean=20_000_000.0,
    )
    print(
        json.dumps(
            {"layer_fraction": layer_fraction.tolist()},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
