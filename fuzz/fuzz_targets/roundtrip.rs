#![no_main]
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(s) = core::str::from_utf8(data) else {
        return;
    };
    let Ok(inchi) = s.parse::<inchi_parser::inchi::InChI>() else {
        return;
    };
    let displayed = inchi.to_string();
    let reparsed = displayed
        .parse::<inchi_parser::inchi::InChI>()
        .unwrap_or_else(|e| {
            panic!("Reparse failed: {displayed}\nOriginal: {s}\nError: {e:?}")
        });
    assert_eq!(
        inchi, reparsed,
        "Round-trip mismatch:\n  original: {s}\n  displayed: {displayed}"
    );
});
