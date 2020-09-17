use std::collections::HashMap;
use std::{env, fs, io};

#[derive(Debug)]
struct Interval(u32, u32, u32, u32);

#[derive(Debug)]
struct IITree {
    a: Vec<Interval>,
}

impl IITree {
    fn new() -> Self {
        IITree { a: Vec::new() }
    }

    fn add(&mut self, start: u32, end: u32, max: u32) {
        self.a.push(Interval(start, end, max, 0))
    }

    fn index(&mut self) {
        self.a.sort_unstable_by_key(|k| k.0);
        let (mut last, mut last_i, mut i) = (0, 1, 0);
        while i < self.a.len() {
            last = self.a[i].1;
            self.a[i].2 = last;
            last_i = i;
            i += 2;
        }
        let mut k = 1;
        while 1 << k <= self.a.len() {
            let (mut i0, step) = (1 << k, 1 << (k + 1));
            while i0 <= self.a.len() {
                let x = 1 << (k - 1);
                let mut max_end = std::cmp::max(self.a[i].1, self.a[i - x].2);
                let e = if i + x <= self.a.len() {
                    self.a[i + x].1
                } else {
                    last
                };
                max_end = std::cmp::max(max_end, e);
                self.a[i] = Interval(self.a[i].3, self.a[i].0, self.a[i].1, max_end);
                i0 += step;
            }
            last_i = if (last_i >> k & 1) != 0 {
                last_i + (1 << (k - 1))
            } else {
                last_i - (1 << (k - 1))
            };

            if last_i <= self.a.len() {
                last = std::cmp::max(last, self.a[last_i].2);
            }
            k += 1;
        }
    }
}

#[derive(Debug)]
struct BedMap<'a> {
    map: HashMap<&'a str, IITree>,
}

impl<'a> BedMap<'a> {
    fn from_str(content: &'a str) -> Self {
        let mut map: HashMap<&str, IITree> = HashMap::with_capacity(24);
        for line in content.lines() {
            let (chrom, start, end) = parse_line(line).expect("Failed to parse BED row.");
            let tree = map.entry(chrom).or_insert(IITree::new());
            tree.add(start, end, end);
        }
        Self { map }
    }

    fn index_trees(&mut self) {
        // sort each grouping
        for (_, tree) in self.map.iter_mut() {
            tree.index();
        }
    }
}

fn parse_line(line: &str) -> Option<(&str, u32, u32)> {
    let mut iter = line.trim().split("\t").take(3);
    let chrom = iter.next()?;
    let start = iter.next()?.parse().ok()?;
    let end = iter.next()?.parse().ok()?;
    Some((chrom, start, end))
}

fn main() -> io::Result<()> {
    let mut args = env::args();

    let path = args.nth(1).expect("Missing input file.");
    let contents = fs::read_to_string(path)?;
    let mut tree = BedMap::from_str(&contents);

    tree.index_trees();
    println!("{:#?}", tree);

    // let path = args.nth(2).expect("Missing other file.");
    // println!("{}", path);

    Ok(())
}

// fn read_bed_by_line() -> io::Result<()> {
//     let path = std::env::args().nth(1).unwrap();
//
//     let f = fs::File::open(path)?;
//     let mut reader = io::BufReader::new(f);
//
//     let mut line = String::new();
//     while let Ok(_) = reader.read_line(&mut line) {
//         let (chrom, start, end) = parse_line(&line).expect("Failed to parse row.");
//         // do something with line
//
//         line.clear();
//     }
//
//     Ok(())
// }
