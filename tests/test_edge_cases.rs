//! Edge case tests for InChI parsing based on known difficult cases
//! from PubChem, RDKit issues, and the InChI specification.

use inchi_parser::inchi::InChI;

fn parse_ok(inchi: &str) {
    assert!(
        inchi.parse::<InChI>().is_ok(),
        "Failed to parse: {inchi}\nError: {:?}",
        inchi.parse::<InChI>().unwrap_err()
    );
}

// ── Proton-only (no formula layer) ──────────────────────────────────────────

#[test]
fn proton() {
    parse_ok("InChI=1S/p+1");
}

#[test]
fn deuteron() {
    parse_ok("InChI=1S/p+1/i/hD");
}

#[test]
fn triton() {
    parse_ok("InChI=1S/p+1/i/hT");
}

#[test]
fn protium_isotope() {
    parse_ok("InChI=1S/p+1/i/hH");
}

#[test]
fn multi_proton() {
    parse_ok("InChI=1S/p+2");
    parse_ok("InChI=1S/p+5");
}

// ── Formula-only (no connections or hydrogen layer) ─────────────────────────

#[test]
fn noble_gases() {
    parse_ok("InChI=1S/He");
    parse_ok("InChI=1S/Ne");
    parse_ok("InChI=1S/Ar");
}

#[test]
fn single_atoms() {
    parse_ok("InChI=1S/H");
    parse_ok("InChI=1S/O");
}

#[test]
fn oxide_anion() {
    parse_ok("InChI=1S/O/q-2");
}

#[test]
fn no_hydrogen_layer() {
    parse_ok("InChI=1S/C2N2/c3-1-2-4");
    parse_ok("InChI=1S/N2/c1-2");
}

#[test]
fn diazene() {
    parse_ok("InChI=1S/H2N2/c1-2/h1-2H");
}

// ── Multi-component with n* repetition ──────────────────────────────────────

#[test]
fn ferrocene() {
    parse_ok("InChI=1S/2C5H5.Fe/c2*1-2-4-5-3-1;/h2*1-5H;");
}

#[test]
fn tetraethyllead() {
    parse_ok("InChI=1S/4C2H5.Pb/c4*1-2;/h4*1H2,2H3;");
}

#[test]
fn cisplatin() {
    parse_ok("InChI=1S/2ClH.2H3N.Pt/h2*1H;2*1H3;/q;;;;+2/p-2");
}

#[test]
fn sodium_acetate() {
    parse_ok("InChI=1S/C2H4O2.Na/c1-2(3)4;/h1H3,(H,3,4);/q;+1/p-1");
}

#[test]
fn sodium_chloride() {
    parse_ok("InChI=1S/ClH.Na/h1H;/q;+1/p-1");
}

#[test]
fn formic_acid_dimer() {
    parse_ok("InChI=1S/2CH2O2/c2*2-1-3/h2*1H,(H,2,3)");
}

#[test]
fn hexafluoroacetone_trihydrate() {
    parse_ok("InChI=1S/C3F6O.3H2O/c4-2(5,6)1(10)3(7,8)9;;;/h;3*1H2");
}

// ── Complex stereochemistry ─────────────────────────────────────────────────

#[test]
fn taxol() {
    parse_ok(
        "InChI=1S/C47H51NO14/c1-25-31(60-43(56)36(52)35(28-16-10-7-11-17-28)48-41(54)29-18-12-8-13-19-29)23-47(57)40(61-42(55)30-20-14-9-15-21-30)38-45(6,32(51)22-33-46(38,24-58-33)62-27(3)50)39(53)37(59-26(2)49)34(25)44(47,4)5/h7-21,31-33,35-38,40,51-52,57H,22-24H2,1-6H3,(H,48,54)/t31-,32-,33+,35-,36+,37+,38-,40-,45+,46-,47+/m0/s1",
    );
}

#[test]
fn morphine() {
    parse_ok(
        "InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1",
    );
}

#[test]
fn borneol_enantiomers() {
    parse_ok("InChI=1S/C10H18O/c1-9(2)7-4-5-10(9,3)8(11)6-7/h7-8,11H,4-6H2,1-3H3/t7-,8+,10+/m0/s1");
    parse_ok("InChI=1S/C10H18O/c1-9(2)7-4-5-10(9,3)8(11)6-7/h7-8,11H,4-6H2,1-3H3/t7-,8+,10+/m1/s1");
}

// ── Mobile hydrogen edge cases ──────────────────────────────────────────────

#[test]
fn urea_all_mobile_h() {
    parse_ok("InChI=1S/CH4N2O/c2-1(3)4/h(H4,2,3,4)");
}

#[test]
fn glycine_mixed_h() {
    parse_ok("InChI=1S/C2H5NO2/c3-1-2(4)5/h1,3H2,(H,4,5)");
}

#[test]
fn fluorosulfonate_mobile_h_with_proton() {
    parse_ok("InChI=1S/FHO3S/c1-5(2,3)4/h(H,2,3,4)/p-1");
}

// ── Isotope layer edge cases ────────────────────────────────────────────────

#[test]
fn carbon_13_methane() {
    parse_ok("InChI=1S/CH4/h1H4/i1+1");
}

#[test]
fn oxygen_18_water() {
    parse_ok("InChI=1S/H2O/h1H2/i1+2");
}

#[test]
fn heavy_water_d2o() {
    parse_ok("InChI=1S/H2O/h1H2/i/hD2");
}

#[test]
fn superheavy_water_t2o() {
    parse_ok("InChI=1S/H2O/h1H2/i/hT2");
}

#[test]
fn semiheavy_water_hdo() {
    parse_ok("InChI=1S/H2O/h1H2/i/hD");
}

#[test]
fn deuterium_gas_combined_shift() {
    parse_ok("InChI=1S/H2/h1H/i1+1D");
}

// ── Isotope before stereo (non-canonical ordering) ──────────────────────────

#[test]
fn isotope_before_double_bond_stereo() {
    parse_ok("InChI=1S/C5H10/c1-4-5(2)3/h4H,1-3H3/i2D3/b5-4-");
}

#[test]
fn isotope_before_tetrahedral_stereo() {
    parse_ok("InChI=1S/C8H17N3/c1-2-3-4-5-6-7-8-10-11-9/h2-8H2,1H3/i8D/t8-/m1/s1");
}

// ── RDKit crashers (syntactically valid) ────────────────────────────────────

#[test]
fn organotin_7_components() {
    parse_ok(
        "InChI=1S/2C18H36N3OP.2CH3.2ClH.Sn/c2*22-23(19-16-10-4-1-5-11-16,20-17-12-6-2-7-13-17)21-18-14-8-3-9-15-18;;;;;/h2*16-18H,1-15H2,(H3,19,20,21,22);2*1H3;2*1H;/q;;2*-1;;;+4/p-2",
    );
}

#[test]
fn organosilicon() {
    parse_ok("InChI=1S/C5H14O3SSi/c1-5-9(10,6-2,7-3)8-4/h5H2,1-4H3");
}

// ── Miscellaneous ───────────────────────────────────────────────────────────

#[test]
fn water() {
    parse_ok("InChI=1S/H2O/h1H2");
}

#[test]
fn methane() {
    parse_ok("InChI=1S/CH4/h1H4");
}

#[test]
fn benzene() {
    parse_ok("InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H");
}

#[test]
fn diborane() {
    // B2H6 has bridging hydrogens numbered as atoms 3-4 in the connection layer.
    parse_ok("InChI=1S/B2H6/c1-3-2-4-1/h1-2H2");
}

#[test]
fn ammonium_nitrate_disconnected() {
    parse_ok("InChI=1S/NO3.H3N/c2-1(3)4;/h;1H3");
}

#[test]
fn incorrect_proton_representation() {
    // InChI=1S/H/q+1 is a valid InChI (single H atom with +1 charge)
    // but NOT the correct way to represent a proton (which is InChI=1S/p+1).
    parse_ok("InChI=1S/H/q+1");
}

// ── Radicals (valid InChIs, no special layer) ───────────────────────────────

#[test]
fn nitric_oxide_radical() {
    parse_ok("InChI=1S/NO/c1-2");
}

#[test]
fn nitrogen_dioxide_radical() {
    parse_ok("InChI=1S/NO2/c2-1-3");
}

#[test]
fn methyl_radical() {
    parse_ok("InChI=1S/CH3/h1H3");
}

// ── Negative proton counts ──────────────────────────────────────────────────

#[test]
fn sodium_hydroxide_p_minus_1() {
    parse_ok("InChI=1S/Na.H2O/h;1H2/q+1;/p-1");
}

#[test]
fn calcium_chloride_p_minus_2() {
    parse_ok("InChI=1S/Ca.2ClH/h;2*1H/q+2;;/p-2");
}

#[test]
fn europium_chloride_p_minus_3() {
    parse_ok("InChI=1S/3ClH.Eu/h3*1H;/q;;;+3/p-3");
}

// ── Multi-charge components with n* in /q ───────────────────────────────────

#[test]
fn potassium_dichromate_11_components() {
    parse_ok("InChI=1S/2Cr.2K.7O/q;;2*+1;;;;;;2*-1");
}

#[test]
fn uranyl_nitrate_actinide() {
    parse_ok("InChI=1S/2NO3.2O.U/c2*2-1(3)4;;;/q2*-1;;;+2");
}

// ── Lanthanides and actinides ───────────────────────────────────────────────

#[test]
fn gadolinium_chloride() {
    parse_ok("InChI=1S/3ClH.Gd/h3*1H;/q;;;+3/p-3");
}

// ── Stereo diastereomers ────────────────────────────────────────────────────

#[test]
fn bromochlorofluoromethane_enantiomers() {
    parse_ok("InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m0/s1");
    parse_ok("InChI=1S/CHBrClF/c2-1(3)4/h1H/t1-/m1/s1");
}

#[test]
fn diastereomers_meso() {
    parse_ok("InChI=1S/C2H2BrClFI/c3-1(6)2(4)5/h1-2H/t1-,2+/m1/s1");
    parse_ok("InChI=1S/C2H2BrClFI/c3-1(6)2(4)5/h1-2H/t1-,2+/m0/s1");
    parse_ok("InChI=1S/C2H2BrClFI/c3-1(6)2(4)5/h1-2H/t1-,2-/m0/s1");
    parse_ok("InChI=1S/C2H2BrClFI/c3-1(6)2(4)5/h1-2H/t1-,2-/m1/s1");
}

// ── Multiple mobile hydrogen groups ─────────────────────────────────────────

#[test]
fn multiple_mobile_h_groups() {
    // Two separate mobile H groups in one /h layer
    parse_ok("InChI=1S/C3H5NO3/c5-2-4-1-3(6)7/h2H,1H2,(H,4,5)(H,6,7)");
}

#[test]
fn large_mobile_h_group() {
    // 6 mobile hydrogens shared among 5 atoms
    parse_ok("InChI=1S/C2H6N4O/c3-1(4)6-2(5)7/h(H6,3,4,5,6,7)");
}

// ── No connection layer (CO2) ───────────────────────────────────────────────

#[test]
fn co2_no_hydrogen_layer() {
    parse_ok("InChI=1S/CO2/c2-1-3");
}

// ── L-glucose (5 stereocenters) ─────────────────────────────────────────────

#[test]
fn l_glucose() {
    parse_ok("InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6+/m1/s1");
}

// ── Silver (single heavy atom, no layers) ───────────────────────────────────

#[test]
fn silver_atom() {
    parse_ok("InChI=1S/Ag");
}
