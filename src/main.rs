use std::{fs::File, io::BufWriter, path::Path, process};

use clap::Parser;
use groan_rs::{
    errors::{GroupError, SimBoxError},
    prelude::*,
};

use std::io::Write;

// Calculate membrane thickness.
#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about,
    long_about = "Calculate a 2D map of membrane thickness."
)]
pub struct Args {
    #[arg(
        short = 's',
        long = "structure",
        help = "Input structure file",
        long_help = "Path to a gro, pdb, or tpr file containing the system structure."
    )]
    structure: String,

    #[arg(
        short = 'f',
        long = "trajectory",
        help = "Input trajectory file",
        long_help = "Path to an xtc file containing the trajectory to analyze."
    )]
    trajectory: String,

    #[arg(
        short = 'o',
        long = "output",
        help = "Path to the output file.",
        long_help = "Path to the output file where the density will be written.",
        default_value = "membrane_thickness.dat"
    )]
    output: String,

    #[arg(
        short = 'n',
        long = "index",
        help = "Input index file",
        long_help = "Path to an ndx file containing groups associated with the system."
    )]
    index: Option<String>,

    #[arg(
        short = 'l',
        long = "lipids",
        help = "Specification of membrane lipids.",
        long_help = "Specify atoms corresponding to membrane lipids.",
        default_value = "@membrane"
    )]
    lipids: String,

    #[arg(
        short = 'p',
        long = "phosphates",
        help = "Specification of atoms identifying lipid headgroups.",
        long_help = "Specify atoms identifying lipid headgroups. Use only one atom per lipid molecule!",
        default_value = "name PO4 P"
    )]
    phosphates: String,

    #[arg(
        short = 'a',
        long = "nan",
        help = "Minimal required number of samples in a grid bin.",
        long_help = "How many phosphates must be detected in a grid bin to calculate membrane thickness for this bin.",
        default_value = "30"
    )]
    nan_limit: usize,

    #[arg(
        long = "xmin",
        help = "Minimum coordinate for the x-dimension of the grid.",
        long_help = "Minimum coordinate for the x-dimension of the grid."
    )]
    xmin: Option<f32>,

    #[arg(
        long = "xmax",
        help = "Maximum coordinate for the x-dimension of the grid.",
        long_help = "Maximum coordinate for the x-dimension of the grid."
    )]
    xmax: Option<f32>,

    #[arg(
        long = "ymin",
        help = "Minimum coordinate for the y-dimension of the grid.",
        long_help = "Minimum coordinate for the y-dimension of the grid."
    )]
    ymin: Option<f32>,

    #[arg(
        long = "ymax",
        help = "Maximum coordinate for the y-dimension of the grid.",
        long_help = "Maximum coordinate for the y-dimension of the grid."
    )]
    ymax: Option<f32>,

    #[arg(
        long = "bin",
        help = "Size of a grid bin in each dimension (in nm).",
        long_help = "Size of a grid bin in each dimension (in nm).",
        default_value_t = 0.1
    )]
    bin_size: f32,
}

/// Print the specified options.
fn print_options(args: &Args, simbox: &SimBox) {
    println!("[STRUCTURE]     {}", args.structure);
    println!("[TRAJECTORY]    {}", args.trajectory);
    println!("[OUTPUT]        {}", args.output);

    if let Some(ndx) = args.index.as_ref() {
        println!("[INDEX]        {}", ndx);
    }

    println!("[LIPIDS]        {}", args.lipids);
    println!("[PHOSPHATES]    {}", args.phosphates);
    println!("[NAN LIMIT]     {}", args.nan_limit);

    println!(
        "[X-RANGE]       {}-{} nm",
        args.xmin.unwrap_or(0.0),
        args.xmax.unwrap_or(simbox.x)
    );
    println!(
        "[Y-RANGE]       {}-{} nm",
        args.ymin.unwrap_or(0.0),
        args.ymax.unwrap_or(simbox.y)
    );

    println!("[BIN SIZE]      {} nm", args.bin_size);
    println!("\n");
}

fn sanity_check_options(args: &Args) -> anyhow::Result<()> {
    if args.nan_limit == 0 {
        anyhow::bail!("NAN limit must be larger than 0, not {}", args.nan_limit);
    }

    if args.xmin > args.xmax {
        anyhow::bail!("Minimum grid x-value cannot be higher than the maximum grid x-value.");
    }

    if args.ymin > args.ymax {
        anyhow::bail!("Minimum grid y-value cannot be higher than the maximum grid y-value.");
    }

    Ok(())
}

fn write_map(
    output_name: impl AsRef<Path>,
    grid_upper: &GridMap<f64, f64, impl Fn(&f64) -> f64>,
    count_upper: &GridMap<usize, usize, impl Fn(&usize) -> usize>,
    grid_lower: &GridMap<f64, f64, impl Fn(&f64) -> f64>,
    count_lower: &GridMap<usize, usize, impl Fn(&usize) -> usize>,
    nan_limit: usize,
    raw_arguments: &Vec<String>,
) -> anyhow::Result<()> {
    let file = File::create(&output_name)?;
    let mut output = BufWriter::new(file);

    writeln!(
        &mut output,
        "# Generated with memthick v{}.",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(&mut output, "# Command line: {}", raw_arguments.join(" "))?;

    writeln!(
        &mut output,
        "# See the average membrane thickness at the end of this file."
    )?;

    writeln!(&mut output, "@ xlabel x-coordinate [nm]")?;
    writeln!(&mut output, "@ ylabel y-coordinate [nm]")?;

    writeln!(&mut output, "@ zlabel membrane thickness [nm]")?;
    writeln!(&mut output, "@ grid --")?;
    writeln!(&mut output, "$ type colorbar")?;
    writeln!(&mut output, "$ colormap rainbow")?;

    let mut average_thickness = Vec::new();
    for (((upper_sum, upper_count), lower_sum), lower_count) in grid_upper
        .extract_raw()
        .zip(count_upper.extract_raw())
        .zip(grid_lower.extract_raw())
        .zip(count_lower.extract_raw())
    {
        let thickness = if *upper_count.2 < nan_limit || *lower_count.2 < nan_limit {
            f64::NAN
        } else {
            let upper_av = upper_sum.2 / (*upper_count.2 as f64);
            let lower_av = lower_sum.2 / (*lower_count.2 as f64);
            upper_av - lower_av
        };

        writeln!(
            &mut output,
            "{:12.6} {:12.6} {:12.4}",
            upper_sum.0, upper_sum.1, thickness
        )?;

        if thickness.is_finite() {
            average_thickness.push(thickness);
        }
    }

    writeln!(
        &mut output,
        "# Average membrane thickness: {:12.4} nm",
        average_thickness.iter().sum::<f64>() / average_thickness.len() as f64
    )?;

    Ok(())
}

fn run() -> anyhow::Result<()> {
    let raw_arguments = std::env::args().collect::<Vec<_>>();

    let args = Args::parse();
    println!("\n>> memthick {} <<\n", env!("CARGO_PKG_VERSION"));
    sanity_check_options(&args)?;

    let mut system = System::from_file(&args.structure).map_err(anyhow::Error::from_boxed)?;
    if let Some(ndx) = &args.index {
        system.read_ndx(ndx)?;
    }

    let simbox = system.get_box().ok_or(SimBoxError::DoesNotExist)?;
    if !simbox.is_orthogonal() || simbox.is_zero() {
        return Err(SimBoxError::NotOrthogonal.into());
    }

    print_options(&args, simbox);

    let xmin = args.xmin.unwrap_or(0.0);
    let xmax = args.xmax.unwrap_or(simbox.x);
    let ymin = args.ymin.unwrap_or(0.0);
    let ymax = args.ymax.unwrap_or(simbox.y);

    match system.group_create("xxxMemthickReservedxxx-Lipids", &args.lipids) {
        Ok(_) | Err(GroupError::AlreadyExistsWarning(_)) => (),
        Err(e) => return Err(e.into()),
    }

    if system
        .group_get_n_atoms("xxxMemthickReservedxxx-Lipids")
        .unwrap()
        == 0
    {
        anyhow::bail!("The query '{}' selects no atoms.", &args.lipids);
    }

    match system.group_create("xxxMemthickReservedxxx-Heads", &args.phosphates) {
        Ok(_) | Err(GroupError::AlreadyExistsWarning(_)) => (),
        Err(e) => return Err(e.into()),
    }

    if system
        .group_get_n_atoms("xxxMemthickReservedxxx-Heads")
        .unwrap()
        == 0
    {
        anyhow::bail!("The query '{}' selects no atoms.", &args.phosphates);
    }

    let mut grid_upper = GridMap::new(
        (xmin, xmax),
        (ymin, ymax),
        (args.bin_size, args.bin_size),
        f64::clone,
    )?;

    let mut grid_lower = GridMap::new(
        (xmin, xmax),
        (ymin, ymax),
        (args.bin_size, args.bin_size),
        f64::clone,
    )?;

    let mut count_upper = GridMap::new(
        (xmin, xmax),
        (ymin, ymax),
        (args.bin_size, args.bin_size),
        usize::clone,
    )?;

    let mut count_lower = GridMap::new(
        (xmin, xmax),
        (ymin, ymax),
        (args.bin_size, args.bin_size),
        usize::clone,
    )?;

    for frame in system
        .group_xtc_iter(&args.trajectory, "xxxMemthickReservedxxx-Lipids")?
        .print_progress(ProgressPrinter::default())
    {
        let frame = frame?;

        let membrane_center = frame
            .group_get_center("xxxMemthickReservedxxx-Lipids")
            .unwrap();

        for head in frame.group_iter("xxxMemthickReservedxxx-Heads").unwrap() {
            let zdist = head
                .distance_from_point(&membrane_center, Dimension::Z, frame.get_box().unwrap())
                .unwrap();

            let position = head.get_position().unwrap();

            let tile_wrapped = if zdist > 0.0 {
                grid_upper.get_mut_at(position.x, position.y)
            } else {
                grid_lower.get_mut_at(position.x, position.y)
            };

            if let Some(tile) = tile_wrapped {
                *tile += zdist as f64;
            }

            let count_wrapped = if zdist > 0.0 {
                count_upper.get_mut_at(position.x, position.y)
            } else {
                count_lower.get_mut_at(position.x, position.y)
            };

            if let Some(count) = count_wrapped {
                *count += 1;
            }
        }
    }

    write_map(
        &args.output,
        &grid_upper,
        &count_upper,
        &grid_lower,
        &count_lower,
        args.nan_limit,
        &raw_arguments,
    )?;

    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{}", e);
        process::exit(1);
    } else {
        process::exit(0);
    }
}
