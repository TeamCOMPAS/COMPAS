from .gpu_utils import xp


def m1_m2_to_chirp_mass(m1, m2):
    return (m1 * m2) ** (3 / 5) / (m1 + m2) ** (1 / 5)


def m1_m2_to_eta(m1, m2):
    return m1 * m2 / (m1 + m2) ** 2


def m1_m2_to_eta_chirp_mass(m1, m2):
    return m1_m2_to_eta(m1, m2), m1_m2_to_chirp_mass(m1, m2)


def chirp_mass_eta_to_total_mass(chirp_mass, eta):
    return chirp_mass / eta ** (3 / 5)


def total_mass_eta_to_m1_m2(total_mass, eta):
    m1 = total_mass * 0.5 * (1. + xp.sqrt(1. - 4 * eta))
    m2 = total_mass - m1
    return m1, m2


def chirp_mass_eta_to_m1_m2(chirp_mass, eta):
    total_mass = chirp_mass_eta_to_total_mass(chirp_mass, eta)
    return total_mass_eta_to_m1_m2(total_mass, eta)
