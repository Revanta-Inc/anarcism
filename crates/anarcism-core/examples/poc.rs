use anarcism_core::{NumberingOptions, number_sequence};

const VH: &str = "EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSAA";
const VL: &str = "DIVMTQSQKFMSTSVGDRVSITCKASQNVGTAVAWYQQKPGQSPKLMIYSASNRYTGVPDRFTGSGSGTDFTLTISNMQSEDLADYFCQQYSSYPLTFGAGTKLELKR";

fn main() {
    for (id, sequence) in [("VH", VH), ("VL", VL)] {
        let result = number_sequence(sequence, &NumberingOptions::default()).unwrap();
        for domain in result.domains {
            println!(
                "{id}\t{:?}\t{}\t{}\t{}\t{:.3}\t{}",
                domain.chain_type,
                domain.species,
                domain.start,
                domain.end,
                domain.bit_score,
                domain.padded_imgt_alignment
            );
            println!(
                "labels\t{}",
                domain
                    .numbering
                    .iter()
                    .map(|residue| format!(
                        "{}{}:{}",
                        residue.position, residue.insertion_code, residue.amino_acid
                    ))
                    .collect::<Vec<_>>()
                    .join(" ")
            );
            println!(
                "alternatives\t{}",
                domain
                    .alternative_hits
                    .iter()
                    .map(|hit| format!("{}:{:.3}", hit.profile, hit.bit_score))
                    .collect::<Vec<_>>()
                    .join(" ")
            );
        }
    }
}
