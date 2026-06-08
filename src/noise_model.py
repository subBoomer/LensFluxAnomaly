import numpy as np

RADIO_SNR = 50.0


class RadioNoise:
    def __init__(self):
        self.snr = RADIO_SNR

    def apply(self, F_true: np.ndarray, rng: np.random.Generator) -> np.ndarray:
        sigma = np.abs(F_true) / self.snr
        return rng.normal(F_true, sigma)
