/// The prefix string for every InChI.
pub const INCHI_PREFIX: &str = "InChI=";

/// Valid single-character prefixes for InChI layers (after the main layer).
///
/// `c` = connections, `h` = H atoms, `q` = charge, `p` = protons,
/// `b` = double bond stereo, `t` = tetrahedral stereo, `m` = stereo type,
/// `s` = stereo bond notation, `i` = isotopes, `f` = fixed H, `r` = reconnected.
pub const KNOWN_LAYER_PREFIXES: &[char] =
    &['c', 'h', 'q', 'p', 'b', 't', 'm', 's', 'i', 'f', 'r'];
