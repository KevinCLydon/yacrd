/*
Copyright (c) 2018 Pierre Marijon <pierre.marijon@inria.fr>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* crates use */
use xz2;
use bzip2;
use flate2;

/* standard use */
use std::io;
use std::fs::File;

#[derive(Debug)]
pub enum CompressionFormat {
    Gzip = 0x1F8B,
    Bzip = 0x425A,
    Lzma = 0xFD377A585A,
    No,
}

pub fn get_input(input_name: &str) -> (Box<io::Read>, CompressionFormat) {
    // choose std::io::stdin or open file
    let raw_input = get_readable(input_name); 

    if input_name == "-" {
        return (Box::new(get_readable(input_name)), CompressionFormat::No);
    }
    // check compression
    let compression = get_compression(raw_input);

    let mut input: Box<io::Read> = get_readable(input_name);
    // return readable and compression status
    match compression {
        CompressionFormat::Gzip => (Box::new(flate2::read::GzDecoder::new(get_readable(input_name))), CompressionFormat::Gzip),
        CompressionFormat::Bzip => (Box::new(bzip2::read::BzDecoder::new(get_readable(input_name))), CompressionFormat::Bzip),
        CompressionFormat::Lzma => (Box::new(xz2::read::XzDecoder::new(get_readable(input_name))), CompressionFormat::Lzma),
        CompressionFormat::No => (Box::new(get_readable(input_name)), CompressionFormat::No),
        
    }
}

fn get_readable(input_name: &str) -> Box<io::Read> {
    match input_name {
        "-" => Box::new(io::stdin()),
        _ => Box::new(File::open(input_name).expect("Can't open input file")),
    }
}

fn get_compression(mut in_stream: Box<io::Read>) -> CompressionFormat {
    let mut buf = vec![0u8; 2];

    in_stream.read_exact(&mut buf);

    match &buf[..] {
        [0x1F, 0x8B] => CompressionFormat::Gzip,
        [0x42, 0x5A] => CompressionFormat::Bzip,
        [0xFD, 0x37] => CompressionFormat::Lzma, // match on 5 value when syntax for subslices in slice are stabilized
        _ => CompressionFormat::No,
    }
}

pub fn get_output(output_name: &str, format: CompressionFormat) -> Box<io::Write> {
    match format {
        CompressionFormat::Gzip => Box::new(flate2::write::GzEncoder::new(get_writable(output_name), flate2::Compression::best())),
        CompressionFormat::Bzip => Box::new(bzip2::write::BzEncoder::new(get_writable(output_name), bzip2::Compression::Best)),
        CompressionFormat::Lzma => Box::new(xz2::write::XzEncoder::new(get_writable(output_name), 9)),
        CompressionFormat::No => Box::new(get_writable(output_name))
    }
}

pub fn choose_compression(input_compression: CompressionFormat, compression_set: bool, compression_value: &str) -> CompressionFormat {
    if !compression_set {
        return input_compression;
    }

    match compression_value {
        "gzip" => CompressionFormat::Gzip,
        "bzip" => CompressionFormat::Bzip,
        "lzma" => CompressionFormat::Lzma,
        _ => CompressionFormat::No,
    }
}

fn get_writable(output_name: &str) -> Box<io::Write> {
    match output_name {
        "-" => Box::new(io::stdout()),
        _ => Box::new(File::create(output_name).expect("Can't open output file")),
    }
}

