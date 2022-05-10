
/* crate use */
use anyhow::{anyhow, bail, Context, Result};
use log::error;

/* local use */
use crate::editor;
use crate::error;
use crate::stack;
use crate::util;

pub fn pluck(
    input_path: &str,
    output_path: &str,
    minimum_bad_region_length: Option<u32>,
    badregions: &mut dyn stack::BadPart,
    not_covered: f64,
    buffer_size: usize,
) -> Result<()> {
    let (input, compression) = util::read_file(input_path, buffer_size)?;
    let output = util::write_file(output_path, compression, buffer_size)?;

    // Default to 50 for minimum_bad_region_length
    let minimum_bad_region_length: u32 = minimum_bad_region_length.unwrap_or(50);

    match util::get_file_type(input_path) {
        Some(util::FileType::Fasta) => fasta(input, output, minimum_bad_region_length, badregions, not_covered)
            .with_context(|| anyhow!("Filename: {}", input_path.to_string()))?,
        Some(util::FileType::Fastq) => fastq(input, output, minimum_bad_region_length, badregions, not_covered)
            .with_context(|| anyhow!("Filename: {}", input_path.to_string()))?,
        Some(util::FileType::Paf) => bail!(error::Error::CantRunOperationOnFile {
            operation: "plucking".to_string(),
            filetype: util::FileType::Paf,
            filename: input_path.to_string()
        }),
        Some(util::FileType::M4) => bail!(error::Error::CantRunOperationOnFile {
            operation: "plucking".to_string(),
            filetype: util::FileType::M4,
            filename: input_path.to_string()
        }),
        Some(util::FileType::Yacrd) => bail!(error::Error::CantRunOperationOnFile {
            operation: "plucking".to_string(),
            filetype: util::FileType::Yacrd,
            filename: input_path.to_string()
        }),
        None | Some(util::FileType::YacrdOverlap) => {
            bail!(error::Error::UnableToDetectFileFormat {
                filename: input_path.to_string()
            })
        }
    };

    Ok(())
}

fn fasta<R, W>(
    input: R,
    output: W,
    minimum_bad_region_length: u32,
    badregions: &mut dyn stack::BadPart,
    not_covered: f64,
) -> Result<()>
    where
        R: std::io::Read,
        W: std::io::Write,
{
    let mut reader = noodles::fasta::Reader::new(std::io::BufReader::new(input));
    let mut writer = noodles::fasta::Writer::new(std::io::BufWriter::new(output));

    for result in reader.records() {
        let record = result.with_context(|| error::Error::ReadingErrorNoFilename {
            format: util::FileType::Fasta,
        })?;

        let (badregion, length) = badregions.get_bad_part(record.name())?;

        let rtype = editor::type_of_read(*length, badregion, not_covered);

        if rtype == editor::ReadType::NotCovered || badregion.is_empty() {
            continue;
        } else {
            for interval in badregion {
                if interval.0 as usize > record.sequence().len()
                    || interval.1 as usize > record.sequence().len()
                {
                    error!("For read {} pluck position is larger than read, it's strange check your data. For this read, this bad region will be ignored.", record.name());
                    break;
                }

                // Skip any intervals shorter than the minimum
                if (interval.1 - interval.0) < minimum_bad_region_length {
                    continue;
                }

                writer
                    .write_record(&noodles::fasta::Record::new(
                        noodles::fasta::record::Definition::new(
                            &format!("{}_{}_{}", record.name(), interval.0, interval.1),
                            None,
                        ),
                        noodles::fasta::record::Sequence::from(
                            record.sequence().as_ref()[(interval.0 as usize)..(interval.1 as usize)]
                                .to_vec(),
                        ),
                    ))
                    .with_context(|| error::Error::WritingErrorNoFilename {
                        format: util::FileType::Fasta,
                    })?;
            }
        }
    }

    Ok(())
}

fn fastq<R, W>(
    input: R,
    output: W,
    minimum_bad_region_length: u32,
    badregions: &mut dyn stack::BadPart,
    not_covered: f64,
) -> Result<()>
    where
        R: std::io::Read,
        W: std::io::Write,
{
    let mut reader = noodles::fastq::Reader::new(std::io::BufReader::new(input));
    let mut writer = noodles::fastq::Writer::new(std::io::BufWriter::new(output));

    for result in reader.records() {
        let record = result.with_context(|| error::Error::ReadingErrorNoFilename {
            format: util::FileType::Fastq,
        })?;

        let (badregion, length) = badregions.get_bad_part(
            std::str::from_utf8(record.name())?
                .split_ascii_whitespace()
                .next()
                .unwrap(),
        )?;

        let rtype = editor::type_of_read(*length, badregion, not_covered);

        if rtype == editor::ReadType::NotCovered || badregion.is_empty() {
            continue;
        } else {
            let mut sequence_description = std::str::from_utf8(record.name())?.splitn(2, ' ');
            let name = sequence_description.next().unwrap();
            let description = sequence_description.next();

            for interval in badregion {
                if interval.0 as usize > record.sequence().len()
                    || interval.1 as usize > record.sequence().len()
                {
                    error!("For read {} pluck position is larger than read, it's strange check your data. For this read, this bad region will be ignored.", name);
                    break;
                }

                // Skip any intervals shorter than the minimum
                if (interval.1 - interval.0) < minimum_bad_region_length {
                    continue;
                }

                writer
                    .write_record(&noodles::fastq::Record::new(
                        match description {
                            Some(desc) => format!("{}_{}_{} {}", name, interval.0, interval.1, desc),
                            None => format!("{}_{}_{}", name, interval.0, interval.1),
                        }
                            .as_bytes(),
                        record.sequence()[(interval.0 as usize)..(interval.1 as usize)].to_vec(),
                        record.quality_scores()[(interval.0 as usize)..(interval.1 as usize)].to_vec(),
                    ))
                    .with_context(|| error::Error::WritingErrorNoFilename {
                        format: util::FileType::Fasta,
                    })?;
            }
        }
    }

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;

    use crate::stack::BadPart;

    use crate::reads2ovl;
    use crate::reads2ovl::Reads2Ovl;

    const FASTA_FILE: &'static [u8] = b">1
ACTGGGGGGACTGGGGGGACTG
>2
ACTG
>3
ACTG
";

    const FASTA_FILE_PLUCKED: &'static [u8] = b">1_0_4
ACTG
>1_9_13
ACTG
>1_18_22
ACTG
";

    #[test]
    fn fasta_keep_begin_end() -> () {
        let mut ovlst = reads2ovl::FullMemory::new(8192);

        ovlst.add_length("1".to_string(), 22);
        ovlst.add_overlap("1".to_string(), (4, 9)).unwrap();
        ovlst.add_overlap("1".to_string(), (13, 18)).unwrap();

        let mut stack = stack::FromOverlap::new(Box::new(ovlst), 0);

        stack.compute_all_bad_part();

        let mut output: Vec<u8> = Vec::new();
        fasta(FASTA_FILE, &mut output, 2, &mut stack, 0.8).unwrap();

        assert_eq!(FASTA_FILE_PLUCKED, &output[..]);
    }

    const FASTA_FILE_PLUCKED2: &'static [u8] = b">1_4_18
GGGGGACTGGGGGG
";

    #[test]
    fn fasta_keep_middle() -> () {
        let mut ovlst = reads2ovl::FullMemory::new(8192);

        ovlst.add_length("1".to_string(), 22);
        ovlst.add_overlap("1".to_string(), (0, 4)).unwrap();
        ovlst.add_overlap("1".to_string(), (18, 22)).unwrap();

        let mut stack = stack::FromOverlap::new(Box::new(ovlst), 0);

        stack.compute_all_bad_part();

        let mut output: Vec<u8> = Vec::new();
        fasta(FASTA_FILE, &mut output, 2, &mut stack, 0.8).unwrap();

        assert_eq!(FASTA_FILE_PLUCKED2, &output[..]);
    }

    const FASTQ_FILE: &'static [u8] = b"@1
ACTGGGGGGACTGGGGGGACTG
+
??????????????????????
@2
ACTG
+
????
@3
ACTG
+
????
";

    const FASTQ_FILE_PLUCKED: &'static [u8] = b"@1_0_4
ACTG
+
????
@1_9_13
ACTG
+
????
@1_18_22
ACTG
+
????
";

    #[test]
    fn fastq_keep_begin_end() {
        let mut ovlst = reads2ovl::FullMemory::new(8192);

        ovlst.add_length("1".to_string(), 22);
        ovlst.add_overlap("1".to_string(), (4, 9)).unwrap();
        ovlst.add_overlap("1".to_string(), (13, 18)).unwrap();

        let mut stack = stack::FromOverlap::new(Box::new(ovlst), 0);

        stack.compute_all_bad_part();

        let mut output: Vec<u8> = Vec::new();
        fastq(FASTQ_FILE, &mut output, 2, &mut stack, 0.8).unwrap();

        assert_eq!(FASTQ_FILE_PLUCKED, &output[..]);
    }

    const FASTQ_FILE_PLUCKED2: &[u8] = b"@1_4_18
GGGGGACTGGGGGG
+
??????????????
";

    #[test]
    fn fastq_keep_middle() {
        let mut ovlst = reads2ovl::FullMemory::new(8192);

        ovlst.add_length("1".to_string(), 22);
        ovlst.add_overlap("1".to_string(), (0, 4)).unwrap();
        ovlst.add_overlap("1".to_string(), (18, 22)).unwrap();

        let mut stack = stack::FromOverlap::new(Box::new(ovlst), 0);

        stack.compute_all_bad_part();

        let mut output: Vec<u8> = Vec::new();
        fastq(FASTQ_FILE, &mut output, 2, &mut stack, 0.8).unwrap();

        assert_eq!(FASTQ_FILE_PLUCKED2, &output[..]);
    }

    const FASTQ_FILE_PLUCKED3: &[u8] = b"";

    #[test]
    fn fastq_keep_none_too_short() {
        let mut ovlst = reads2ovl::FullMemory::new(8192);

        ovlst.add_length("1".to_string(), 22);
        ovlst.add_overlap("1".to_string(), (0, 4)).unwrap();
        ovlst.add_overlap("1".to_string(), (18, 22)).unwrap();

        let mut stack = stack::FromOverlap::new(Box::new(ovlst), 0);

        stack.compute_all_bad_part();

        let mut output: Vec<u8> = Vec::new();
        fastq(FASTQ_FILE, &mut output, 50, &mut stack, 0.8).unwrap();

        assert_eq!(FASTQ_FILE_PLUCKED3, &output[..]);
    }
}