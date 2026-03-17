//! Display implementations for InChI and its sublayers.

use alloc::{string::String, vec::Vec};
use core::fmt::{self, Display, Formatter, Write};

use elements_rs::isotopes::HydrogenIsotope;

use crate::{
    inchi::{
        InChI,
        charge_layer::ChargeSubLayer,
        isotope_layer::{IsotopeAtom, IsotopeComponent, IsotopeHydrogen, IsotopeLayer},
        main_layer::{HydrogenComponent, HydrogensSubLayer, MainLayer},
        proton_layer::ProtonSublayer,
        stereochemistry_layer::{
            AlleneSublayer, DoubleBondStereo, DoubleBondSublayer,
            StereoChemistryInformationSublayer, StereoParity, StereochemistryLayer,
            TetrahedralComponent, TetrahedralStereo, TetrahedralSublayer,
        },
    },
    version::Version,
};

// ---------------------------------------------------------------------------
// InChI<V>
// ---------------------------------------------------------------------------

impl<V: Version> Display for InChI<V> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}/", crate::constants::INCHI_PREFIX, V::VERSION_PREFIX)?;

        if let Some(main) = &self.main_layer {
            Display::fmt(main, f)?;
        }
        if let Some(charge) = &self.charge {
            Display::fmt(charge, f)?;
        }
        if let Some(proton) = &self.proton
            && proton.proton_count != 0
        {
            if self.main_layer.is_some() {
                Display::fmt(proton, f)?;
            } else {
                // Proton-only: emit `p±N` directly (the `/` was already
                // written as the version separator above).
                f.write_char('p')?;
                if proton.proton_count >= 0 {
                    write!(f, "+{}", proton.proton_count)?;
                } else {
                    write!(f, "{}", proton.proton_count)?;
                }
            }
        }
        if let Some(stereo) = &self.stereochemistry {
            Display::fmt(stereo, f)?;
        }
        if let Some(isotope) = &self.isotope {
            Display::fmt(isotope, f)?;
        }
        if let Some(iso_stereo) = &self.isotope_stereochemistry {
            Display::fmt(iso_stereo, f)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// MainLayer
// ---------------------------------------------------------------------------

impl Display for MainLayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        Display::fmt(&self.chemical_formula, f)?;
        if let Some(connections) = &self.atom_connections {
            write!(f, "/c")?;
            let segments: Vec<String> =
                connections.iter().map(|comp| comp.connection_string.clone()).collect();
            fmt_components_with_repetition(f, &segments, ';')?;
        }
        if let Some(hydrogens) = &self.hydrogens {
            write!(f, "/h")?;
            Display::fmt(hydrogens, f)?;
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Hydrogen layer
// ---------------------------------------------------------------------------

impl Display for HydrogensSubLayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let segments: Vec<String> = self
            .components
            .iter()
            .map(|comp| {
                let mut s = String::new();
                fmt_hydrogen_component(&mut s, comp).unwrap();
                s
            })
            .collect();
        fmt_components_with_repetition(f, &segments, ';')
    }
}

/// Format a single hydrogen component.
fn fmt_hydrogen_component(f: &mut impl Write, comp: &HydrogenComponent) -> fmt::Result {
    let mut first = true;

    // 1) Fixed hydrogens: group atoms by H-count, sort ascending.
    let mut groups: Vec<(u8, Vec<u16>)> = Vec::new();
    for (idx, &count) in comp.fixed_h.iter().enumerate() {
        if count == 0 {
            continue;
        }
        let idx = u16::try_from(idx).expect("atom index fits in u16");
        if let Some(g) = groups.iter_mut().find(|(c, _)| *c == count) {
            g.1.push(idx);
        } else {
            groups.push((count, alloc::vec![idx]));
        }
    }
    groups.sort_by_key(|(count, _)| *count);

    for (count, atoms) in &groups {
        // atoms within each group are already ascending (iterated in order)
        // Merge consecutive atoms into ranges.
        let ranges = merge_ranges(atoms);
        for range in &ranges {
            if !first {
                f.write_char(',')?;
            }
            first = false;
            write!(f, "{}", range.start + 1)?;
            if range.start != range.end {
                write!(f, "-{}", range.end + 1)?;
            }
        }
        f.write_char('H')?;
        if *count > 1 {
            write!(f, "{count}")?;
        }
    }

    // 2) Mobile hydrogen groups — no comma between consecutive groups, but a comma
    //    separates fixed-H entries from the first mobile group.
    let mut prev_was_mobile = false;
    for mobile in &comp.mobile_groups {
        if !first && !prev_was_mobile {
            f.write_char(',')?;
        }
        first = false;
        prev_was_mobile = true;
        f.write_char('(')?;
        f.write_char('H')?;
        if mobile.count > 1 {
            write!(f, "{}", mobile.count)?;
        }
        if mobile.negative_count > 0 {
            f.write_char('-')?;
            if mobile.negative_count > 1 {
                write!(f, "{}", mobile.negative_count)?;
            }
        }
        for atom in &mobile.atoms {
            write!(f, ",{}", atom + 1)?;
        }
        f.write_char(')')?;
    }

    Ok(())
}

/// Merge sorted atom indices into inclusive ranges of consecutive values.
struct Range {
    start: u16,
    end: u16,
}

fn merge_ranges(sorted_atoms: &[u16]) -> Vec<Range> {
    let mut ranges: Vec<Range> = Vec::new();
    for &atom in sorted_atoms {
        if let Some(last) = ranges.last_mut()
            && atom == last.end + 1
        {
            last.end = atom;
            continue;
        }
        ranges.push(Range { start: atom, end: atom });
    }
    ranges
}

// ---------------------------------------------------------------------------
// Repetition helper
// ---------------------------------------------------------------------------

/// Formats segments joined by `separator`, compressing consecutive identical
/// non-empty segments with an `n*` repetition prefix.
fn fmt_components_with_repetition(
    f: &mut Formatter<'_>,
    segments: &[String],
    separator: char,
) -> fmt::Result {
    let mut i = 0;
    while i < segments.len() {
        if i > 0 {
            f.write_char(separator)?;
        }
        if segments[i].is_empty() {
            i += 1;
        } else {
            let mut run_len = 1;
            while i + run_len < segments.len() && segments[i + run_len] == segments[i] {
                run_len += 1;
            }
            if run_len > 1 {
                write!(f, "{run_len}*")?;
            }
            f.write_str(&segments[i])?;
            i += run_len;
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Charge layer
// ---------------------------------------------------------------------------

impl Display for ChargeSubLayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("/q")?;
        let segments: Vec<String> = self
            .charges
            .iter()
            .map(|&charge| {
                match charge.cmp(&0) {
                    core::cmp::Ordering::Greater => alloc::format!("+{charge}"),
                    core::cmp::Ordering::Less => alloc::format!("{charge}"),
                    core::cmp::Ordering::Equal => String::new(),
                }
            })
            .collect();
        fmt_components_with_repetition(f, &segments, ';')
    }
}

// ---------------------------------------------------------------------------
// Proton layer
// ---------------------------------------------------------------------------

impl Display for ProtonSublayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "/p")?;
        if self.proton_count >= 0 {
            write!(f, "+{}", self.proton_count)
        } else {
            write!(f, "{}", self.proton_count)
        }
    }
}

// ---------------------------------------------------------------------------
// Stereochemistry layer
// ---------------------------------------------------------------------------

impl Display for StereochemistryLayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        if let Some(db) = &self.double_bond {
            Display::fmt(db, f)?;
        }
        if let Some(tet) = &self.tetrahedral {
            Display::fmt(tet, f)?;
        }
        if let Some(allene) = &self.allene {
            Display::fmt(allene, f)?;
        }
        if let Some(info) = &self.stereo_info {
            Display::fmt(info, f)?;
        }
        Ok(())
    }
}

fn fmt_parity(f: &mut Formatter<'_>, parity: StereoParity) -> fmt::Result {
    match parity {
        StereoParity::Plus => f.write_char('+'),
        StereoParity::Minus => f.write_char('-'),
        StereoParity::Unknown => f.write_char('?'),
    }
}

impl Display for DoubleBondSublayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("/b")?;
        let segments: Vec<String> = self
            .components
            .iter()
            .map(|comp| {
                let mut s = String::new();
                for (j, bond) in comp.iter().enumerate() {
                    if j > 0 {
                        s.push(',');
                    }
                    write!(s, "{bond}").unwrap();
                }
                s
            })
            .collect();
        fmt_components_with_repetition(f, &segments, ';')
    }
}

impl Display for DoubleBondStereo {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}-{}", self.atom1 + 1, self.atom2 + 1)?;
        fmt_parity(f, self.parity)
    }
}

impl Display for TetrahedralSublayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("/t")?;
        let has_abbreviation =
            self.components.iter().any(|c| matches!(c, TetrahedralComponent::SameAsMainLayer));
        let segments: Vec<String> = self
            .components
            .iter()
            .map(|comp| {
                match comp {
                    TetrahedralComponent::SameAsMainLayer => String::from("m"),
                    TetrahedralComponent::Explicit(centers) => {
                        let mut s = String::new();
                        for (j, center) in centers.iter().enumerate() {
                            if j > 0 {
                                s.push(',');
                            }
                            write!(s, "{center}").unwrap();
                        }
                        s
                    }
                }
            })
            .collect();
        if has_abbreviation {
            // Abbreviations must not be compressed with n* repetition
            for (i, seg) in segments.iter().enumerate() {
                if i > 0 {
                    f.write_char(';')?;
                }
                f.write_str(seg)?;
            }
            Ok(())
        } else {
            fmt_components_with_repetition(f, &segments, ';')
        }
    }
}

impl Display for TetrahedralStereo {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.atom + 1)?;
        fmt_parity(f, self.parity)
    }
}

impl Display for AlleneSublayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("/m")?;
        for (i, group) in self.values.iter().enumerate() {
            if i > 0 {
                f.write_char('.')?;
            }
            for v in group {
                write!(f, "{v}")?;
            }
        }
        Ok(())
    }
}

impl Display for StereoChemistryInformationSublayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "/s{}", self.value)
    }
}

// ---------------------------------------------------------------------------
// Isotope layer
// ---------------------------------------------------------------------------

impl Display for IsotopeLayer {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str("/i")?;
        let segments: Vec<String> = self
            .components
            .iter()
            .map(|comp| {
                let mut s = String::new();
                fmt_isotope_component_atoms(&mut s, comp).unwrap();
                s
            })
            .collect();
        fmt_components_with_repetition(f, &segments, ';')?;
        // Hydrogen isotope sublayer: if any component has hydrogens, emit
        // /h followed by the isotope specs. The /h sublayer is shared across
        // components (component 0's hydrogens are used).
        if let Some(comp) = self.components.first()
            && !comp.hydrogens.is_empty()
        {
            f.write_str("/h")?;
            fmt_isotope_hydrogens(f, &comp.hydrogens)?;
        }
        Ok(())
    }
}

/// Format the atom-level isotope specs for one component.
fn fmt_isotope_component_atoms(f: &mut impl Write, comp: &IsotopeComponent) -> fmt::Result {
    for (i, atom) in comp.atoms.iter().enumerate() {
        if i > 0 {
            f.write_char(',')?;
        }
        fmt_isotope_atom(f, atom)?;
    }
    Ok(())
}

/// Format a single atom isotope specification.
fn fmt_isotope_atom(f: &mut impl Write, atom: &IsotopeAtom) -> fmt::Result {
    write!(f, "{}", atom.atom_index + 1)?;
    if let Some(shift) = atom.mass_shift {
        if shift >= 0 {
            write!(f, "+{shift}")?;
        } else {
            write!(f, "{shift}")?;
        }
    }
    for &h in &atom.hydrogen_isotopes {
        fmt_isotope_hydrogen_inline(f, h)?;
    }
    Ok(())
}

/// Format a hydrogen isotope inline (D/T/H with count).
fn fmt_isotope_hydrogen_inline(f: &mut impl Write, h: IsotopeHydrogen) -> fmt::Result {
    match h.isotope {
        HydrogenIsotope::D => f.write_char('D')?,
        HydrogenIsotope::T => f.write_char('T')?,
        HydrogenIsotope::H1 => f.write_char('H')?,
        _ => unreachable!("InChI only uses D, T, H1 hydrogen isotopes"),
    }
    if h.count > 1 {
        write!(f, "{}", h.count)?;
    }
    Ok(())
}

/// Format the hydrogen isotope sublayer.
fn fmt_isotope_hydrogens(f: &mut impl Write, hydrogens: &[IsotopeHydrogen]) -> fmt::Result {
    for &h in hydrogens {
        fmt_isotope_hydrogen_inline(f, h)?;
    }
    Ok(())
}
