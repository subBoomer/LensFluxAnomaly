SCHEMA = {
    'simulation': {
        'n_systems': int,
        'checkpoint_interval': int,
        'seed': int,
    },
    'population': {
        'z_l_source': str,
        'z_s_min': (int, float),
        'z_s_max': (int, float),
        'sigma_v_mu': (int, float),
        'sigma_v_sigma': (int, float),
        'q_mean': (int, float),
        'q_sigma': (int, float),
        'q_min': (int, float),
        'q_max': (int, float),
        'gamma_ext_scale': (int, float),
    },
    'substructure': {
        'f_sub_dmo': (int, float),
        'f_sub_reference': str,
        'mass_slope': (int, float),
        'mass_min': (int, float),
        'mass_max': (int, float),
        'lambda_dmo': (int, float),
        'concentration_model': str,
    },
    'selection': {
        'survey': str,
        'beta_magnification': (int, float),
        'flux_limit_mjy': (int, float),
        'beam_arcsec': (int, float),
        'quad_detect_prob': (int, float),
    },
    'noise': {
        'radio_snr': (int, float),
    },
    'comparison': {
        'n_bootstrap': int,
    },
}

REQUIRED_META = ['project', 'version', 'cosmology']
COSMOLOGY_KEYS = ['h0', 'omega_m', 'omega_lambda', 'reference']


def validate(config: dict) -> None:
    for key in REQUIRED_META:
        if key not in config:
            raise ValueError(f"config.yaml: missing required top-level key '{key}'")

    cosmo = config.get('cosmology', {})
    for key in COSMOLOGY_KEYS:
        if key not in cosmo:
            raise ValueError(f"config.yaml: cosmology missing required key '{key}'")

    for section, keys in SCHEMA.items():
        if section not in config:
            raise ValueError(f"config.yaml: missing required section '{section}'")
        cfg_sec = config[section]
        for key, expected_type in keys.items():
            if key not in cfg_sec:
                raise ValueError(f"config.yaml: '{section}' missing key '{key}'")
            val = cfg_sec[key]
            if not isinstance(val, expected_type):
                raise ValueError(
                    f"config.yaml: '{section}.{key}' expected {expected_type}, got {type(val).__name__}"
                )
