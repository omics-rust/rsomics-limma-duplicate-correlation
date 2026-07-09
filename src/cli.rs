use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_limma_duplicate_correlation::{Options, run, write_results};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-limma-duplicate-correlation", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    /// log-expression matrix TSV: header = sample ids, col 1 = gene ids.
    pub expr: PathBuf,
    /// Design matrix TSV: header = coefficient names, col 1 = sample ids.
    #[arg(long)]
    design: PathBuf,
    /// Block factor TSV: col 1 = sample ids, col 2 = block label (one row per sample).
    #[arg(long)]
    block: PathBuf,
    /// Result TSV destination; "-" is stdout.
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let opts = Options {
            expr: &self.expr,
            design: &self.design,
            block: &self.block,
            threads: self.common.thread_count(),
        };
        let res = run(&opts)?;

        let mut out: Box<dyn std::io::Write> = if self.output == "-" && self.common.json {
            Box::new(std::io::sink())
        } else if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        write_results(&res, &mut out)?;

        if !self.common.quiet {
            eprintln!(
                "{} genes ({} used), consensus.correlation = {:.6}",
                res.genes.len(),
                res.n_used,
                res.consensus
            );
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "REML consensus intra-block correlation for a log-expression matrix.",
    origin: Some(Origin {
        upstream: "limma duplicateCorrelation",
        upstream_license: "GPL (>=2)",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/bti270"),
    }),
    usage_lines: &["<expr.tsv> --design <design.tsv> --block <blocks.tsv> [-o result.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "design",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Design matrix TSV (header = coefficient names, col 1 = sample ids).",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "block",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Block factor TSV (col 1 = sample ids, col 2 = block label).",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("String"),
                required: false,
                default: Some("-"),
                description: "Result TSV destination; \"-\" is stdout.",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Technical-replicate correlation across blocks",
        command: "rsomics-limma-duplicate-correlation E.tsv --design design.tsv --block blocks.tsv -o cor.tsv",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
