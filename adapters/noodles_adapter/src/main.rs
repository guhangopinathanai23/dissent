//! noodles adapter -- voice #3, an independent pure-Rust VCF parser.
//!
//! Emits the SAME JSON summary contract as the Python adapters so the runner
//! treats every parser identically. Each record is:
//!   [chrom, pos, ref, [alts...], {info}]
//!
//! We use RecordBuf (the materialized, typed record) so INFO comes back typed,
//! the same way pysam/vcfpy hand back typed INFO. INFO Float is f32 in noodles
//! (as in htslib); we widen it to f64 on the way out so it round-trips to the
//! same JSON as pysam's float32-derived value (and the real f32-vs-f64 split is
//! against vcfpy, which uses f64).
//!
//! Canonicalization of missing tokens / unordered INFO / etc. is done centrally
//! in adapters/common.py -- this adapter only faithfully extracts.

use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{self, BufReader, Write};
use std::process::ExitCode;

use noodles::vcf;
use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
use noodles::vcf::variant::record_buf::info::field::value::Array as InfoArray;
use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
use noodles::vcf::variant::record_buf::samples::sample::value::Array as SampleArray;
use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;
use noodles::vcf::variant::record_buf::samples::sample::Value as SampleValue;
use serde_json::{json, Map, Value};

const TOOL: &str = "noodles (rust)";

fn opt_num<T: Into<f64> + Copy>(o: &Option<T>, float: bool) -> Value {
    match o {
        Some(x) => {
            let f: f64 = (*x).into();
            if float { json!(f) } else { json!(f as i64) }
        }
        None => Value::Null,
    }
}

fn info_array_to_json(a: &InfoArray) -> Value {
    match a {
        InfoArray::Integer(xs) => Value::Array(xs.iter().map(|o| opt_num(o, false)).collect()),
        InfoArray::Float(xs) => Value::Array(xs.iter().map(|o| opt_num(o, true)).collect()),
        InfoArray::Character(xs) => Value::Array(
            xs.iter().map(|o| match o { Some(c) => json!(c.to_string()), None => Value::Null }).collect(),
        ),
        InfoArray::String(xs) => Value::Array(
            xs.iter().map(|o| match o { Some(s) => json!(s), None => Value::Null }).collect(),
        ),
    }
}

fn info_value_to_json(v: &InfoValue) -> Value {
    match v {
        InfoValue::Integer(i) => json!(i),
        InfoValue::Float(f) => json!(*f as f64), // widen f32 -> f64 to match htslib/pysam
        InfoValue::Flag => json!(true),
        InfoValue::Character(c) => json!(c.to_string()),
        InfoValue::String(s) => json!(s),
        InfoValue::Array(a) => info_array_to_json(a),
    }
}

fn genotype_to_string(g: &Genotype) -> String {
    // LITERAL reconstruction of the VCF GT string: walk alleles in order, '.'
    // for a missing allele, '|'/'/' per each allele's own phasing. No sorting,
    // no phasing convention imposed -- equality lives in the central canonicalizer.
    let mut s = String::new();
    for (i, allele) in g.as_ref().iter().enumerate() {
        if i > 0 {
            s.push(match allele.phasing() {
                Phasing::Phased => '|',
                Phasing::Unphased => '/',
            });
        }
        match allele.position() {
            Some(p) => s.push_str(&p.to_string()),
            None => s.push('.'),
        }
    }
    s
}

fn sample_array_to_json(a: &SampleArray) -> Value {
    match a {
        SampleArray::Integer(xs) => Value::Array(xs.iter().map(|o| opt_num(o, false)).collect()),
        SampleArray::Float(xs) => Value::Array(xs.iter().map(|o| opt_num(o, true)).collect()),
        SampleArray::Character(xs) => Value::Array(
            xs.iter().map(|o| match o { Some(c) => json!(c.to_string()), None => Value::Null }).collect(),
        ),
        SampleArray::String(xs) => Value::Array(
            xs.iter().map(|o| match o { Some(s) => json!(s), None => Value::Null }).collect(),
        ),
    }
}

fn sample_value_to_json(v: &SampleValue) -> Value {
    match v {
        SampleValue::Integer(i) => json!(i),
        SampleValue::Float(f) => json!(*f as f64),
        SampleValue::Character(c) => json!(c.to_string()),
        SampleValue::String(s) => json!(s),
        SampleValue::Genotype(g) => json!(genotype_to_string(g)),
        SampleValue::Array(a) => sample_array_to_json(a),
    }
}

fn parse_records(path: &str) -> io::Result<Vec<Value>> {
    let file = File::open(path)?;
    let mut reader = vcf::io::Reader::new(BufReader::new(file));
    let header = reader.read_header()?;

    // Only compare INFO fields DECLARED in this file's header (see the Python
    // adapters for the rationale); undeclared fields are excluded.
    let declared: HashSet<String> = header.infos().keys().cloned().collect();

    let mut out: Vec<Value> = Vec::new();
    for result in reader.record_bufs(&header) {
        let record = result?;

        let chrom = record.reference_sequence_name().to_string();
        let pos: u64 = match record.variant_start() {
            Some(p) => p.get() as u64,
            None => 0,
        };
        let reference = record.reference_bases().to_string();
        let alts: Vec<String> = record.alternate_bases().as_ref().to_vec();

        let mut info = Map::new();
        for (key, opt_value) in record.info().as_ref() {
            if !declared.contains(key) {
                continue;
            }
            let jv = match opt_value {
                Some(v) => info_value_to_json(v),
                None => Value::Null,
            };
            info.insert(key.clone(), jv);
        }

        let sample_keys: Vec<String> = record.samples().keys().as_ref().iter().cloned().collect();
        let mut samples: Vec<Value> = Vec::new();
        for sample in record.samples().values() {
            let mut smap = Map::new();
            for (key, opt_value) in sample_keys.iter().zip(sample.values()) {
                let jv = match opt_value {
                    Some(v) => sample_value_to_json(v),
                    None => Value::Null,
                };
                smap.insert(key.clone(), jv);
            }
            samples.push(Value::Object(smap));
        }

        out.push(json!([chrom, pos, reference, alts, Value::Object(info), samples]));
    }
    Ok(out)
}

fn main() -> ExitCode {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: noodles_adapter <vcf_path>");
        return ExitCode::from(2);
    }

    let summary = match parse_records(&args[1]) {
        Ok(records) => json!({
            "tool": TOOL,
            "kind": "ok",
            "num_records": records.len(),
            "records": records,
            "error": Value::Null,
            "notes": Vec::<String>::new(),
        }),
        Err(e) => json!({
            "tool": TOOL,
            "kind": "exception",
            "num_records": 0,
            "records": Vec::<Value>::new(),
            "error": format!("{:?}: {}", e.kind(), e),
        }),
    };

    let mut stdout = io::stdout();
    let _ = stdout.write_all(summary.to_string().as_bytes());
    let _ = stdout.flush();
    ExitCode::SUCCESS
}
