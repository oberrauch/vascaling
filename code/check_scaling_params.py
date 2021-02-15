"""The scaling relations work with SI units (i.e. volume im m3 and area in m2)
or with km3 and km2. While the scaling exponent doesn't change, the scaling
constants are affected, due to the non-linearity of the equation. """

from numpy.testing import assert_allclose


def convert_c_a(c_a, to_SI=True, gamma=1.375):
    return c_a * 10 ** ((1 if to_SI else -1) * 3 * (3 - 2 * gamma))


def convert_c_l(c_l, to_SI=True, q=2.2):
    return c_l * 10 ** ((1 if to_SI else -1) * 3 * (3 - q))


def vascaling(area, c_a, gamma=1.375):
    return c_a * area ** gamma


def compare_volume_area_scaling(c_a_m2, c_a_km2, area=8e6, gamma=1.375):
    v_m3 = vascaling(area, c_a_m2, gamma)
    v_km3 = vascaling(area * 1e-6, c_a_km2, gamma)
    print(f'Volume: {v_m3*1e-9:.3f}x10^9 m3 vs. {v_km3:.3f} km3')
    assert_allclose(v_m3, v_km3*1e9, rtol=1e-4)


def cycle_c_a(c_a, to_SI=True, gamma=1.375):
    print('Converting volume/area scaling constant\n-----------------------')
    unit_in = 'km^(3-2*gamma)' if to_SI else 'm^(3-2*gamma)'
    unit_out = 'm^(3-2*gamma)' if to_SI else 'km^(3-2*gamma)'
    # convert
    c_a_conv = convert_c_a(c_a, to_SI, gamma)
    print(f'From {c_a:.4f} {unit_in} to {c_a_conv:.4f} {unit_out}')
    # convert back
    c_a_back = convert_c_a(c_a_conv, not to_SI, gamma)
    print(f'From {c_a_conv:.4f} {unit_out} back to {c_a_back:.4f} {unit_in}')
    assert_allclose(c_a, c_a_back)
    c_a_m2 = c_a_conv if to_SI else c_a
    c_a_km2 = c_a if to_SI else c_a_conv
    compare_volume_area_scaling(c_a_m2, c_a_km2, gamma=gamma)


def vlscaling(length, c_l, q=2.2):
    return c_l * length ** q


def compare_volume_length_scaling(c_l_m, c_l_km, length=7.1e3, q=2.2):
    v_m3 = vlscaling(length, c_l_m, q)
    v_km3 = vlscaling(length * 1e-3, c_l_km, q)
    print(f'Volume: {v_m3*1e-9:.3f}x10^9 m3 vs. {v_km3:.3f} km3')
    assert_allclose(v_m3, v_km3 * 1e9, rtol=1e-4)


def cycle_c_l(c_l, to_SI=True, q=2.2):
    print('Converting volume/length scaling constant\n-----------------------')
    unit_in = 'km^(3-q)' if to_SI else 'm^(3-q)'
    unit_out = 'm^(3-q)' if to_SI else 'km^(3-q)'
    # convert
    c_l_conv = convert_c_l(c_l, to_SI, q)
    print(f'From {c_l:.4f} {unit_in} to {c_l_conv:.4f} {unit_out}')
    # convert back
    c_l_back = convert_c_l(c_l_conv, not to_SI, q)
    print(f'From {c_l_conv:.4f} {unit_out} back to {c_l_back:.4f} {unit_in}')
    assert_allclose(c_l, c_l_back)
    c_l_m2 = c_l_conv if to_SI else c_l
    c_l_km2 = c_l if to_SI else c_l_conv
    compare_volume_length_scaling(c_l_m2, c_l_km2, q=q)


if __name__ == '__main__':
    cycle_c_a(0.0538, gamma=1.25)
    print()
    cycle_c_l(0.2252, q=2.5)
