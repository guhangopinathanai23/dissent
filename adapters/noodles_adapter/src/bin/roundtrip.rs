//! Round-trip probe (investigation tool, not part of the fuzzer).
//!
//! Reads a VCF with noodles and re-emits each record through noodles' OWN
//! writer, so we can see exactly what bytes noodles would WRITE back. Used to
//! answer: for an empty/`.` allele, do pysam and noodles actually emit
//! different data, or only hold it differently in memory?

use std::env;
use std::fs::File;
use std::io::{self, BufReader, Write};

use noodles::vcf;
use noodles::vcf::variant::io::Write as _;

fn main() -> io::Result<()> {
    let path = env::args().nth(1).expect("usage: roundtrip <vcf>");
    let file = File::open(&path)?;
    let mut reader = vcf::io::Reader::new(BufReader::new(file));
    let header = reader.read_header()?;

    let mut writer = vcf::io::Writer::new(Vec::new());
    for result in reader.records() {
        let record = result?;
        writer.write_variant_record(&header, &record)?;
    }
    io::stdout().write_all(writer.get_ref())?;
    Ok(())
}
