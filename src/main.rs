use std::collections::HashMap;
use std::io::prelude::*;
use std::{env, fs, io};

#[derive(Debug, Clone, Copy)]
struct Interval {
    st: u32,
    en: u32,
    max: u32,
    data: u32,
}

#[derive(Debug)]
struct IITree {
    a: Vec<Interval>,
    stack: Vec<(usize, usize, usize)>,
}

impl IITree {
    fn new() -> Self {
        IITree {
            a: Vec::new(),
            stack: Vec::with_capacity(64),
        }
    }

    fn add(&mut self, st: u32, en: u32, max: u32, data: u32) {
        self.a.push(Interval { st, en, max, data })
    }

    fn index(&mut self) -> usize {
        self.a.sort_unstable_by_key(|k| k.st);

        let (mut last, mut last_i, mut i) = (0, 1, 0);
        while i < self.a.len() {
            last = self.a[i].en;
            self.a[i].max = self.a[i].en;
            last_i = i;
            i += 2;
        }

        let mut k = 1;
        while 1 << k <= self.a.len() {
            let step = 1 << (k + 1);
            let mut i = (1 << k) - 1;
            while i < self.a.len() {
                let x = 1 << (k - 1);
                let max = std::cmp::max(self.a[i].en, self.a[i - x].max);
                let e = if i + x < self.a.len() {
                    self.a[i + x].max
                } else {
                    last
                };
                let max = std::cmp::max(e, max);
                self.a[i].max = max;

                i += step;
            }

            last_i = if (last_i >> k & 1) != 0 {
                last_i - (1 << (k - 1))
            } else {
                last_i + (1 << (k - 1))
            };

            if last_i < self.a.len() {
                last = std::cmp::max(last, self.a[last_i].max);
            }

            k += 1;
        }

        k - 1
    }

    fn overlap(&mut self, st: u32, en: u32, b: &mut Vec<Interval>) {
        let mut h = 0;
        while 1 << h <= self.a.len() {
            h += 1;
        }
        h -= 1;
        self.stack.push(((1 << h) - 1, h, 0));
        while let Some((x, h, w)) = self.stack.pop() {
            if h <= 3 {
                let mut i = x >> h << h;
                let mut i1 = i + (1 << (h + 1)) - 1;
                if i1 >= self.a.len() {
                    i1 = self.a.len();
                }
                while i < i1 && self.a[i].st < en {
                    if st < self.a[i].en {
                        b.push(self.a[i]);
                    }
                    i += 1;
                }
            } else if w == 0 {
                self.stack.push((x, h, 1));
                let y = x - (1 << (h - 1));
                if y >= self.a.len() || self.a[y].max > st {
                    self.stack.push((y, h - 1, 0));
                }
            } else if x < self.a.len() && self.a[x].st < en {
                if st < self.a[x].en {
                    b.push(self.a[x]);
                }
                self.stack.push((x + (1 << (h - 1)), h - 1, 0));
            }
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
    let path = env::args().nth(1).expect("Missing input file.");
    let content = fs::read_to_string(path)?;
    let stdout = io::stdout();
    let mut stdout = stdout.lock();

    let mut map: HashMap<&str, IITree> = HashMap::with_capacity(24);
    for (i, line) in content.lines().enumerate() {
        let (chrom, start, end) = parse_line(line).expect("Failed to parse BED row.");
        let tree = map.entry(chrom).or_insert(IITree::new());
        tree.add(start, end, end, i as u32);
    }

    for (_, tree) in map.iter_mut() {
        tree.index();
    }

    let path = env::args().nth(2).expect("Missing file 2.");
    let f = fs::File::open(path)?;
    let mut reader = io::BufReader::new(f);
    let mut line = String::new();

    let mut b: Vec<Interval> = Vec::with_capacity(200);

    reader
        .read_line(&mut line)
        .expect("Failed to parse record.");
    while !line.is_empty() {
        let (chrom, st0, en0) = parse_line(&line).expect("failed to parse row.");
        let (start, end, len, cov) = match map.get_mut(chrom) {
            Some(tree) => {
                tree.overlap(st0, en0, &mut b);
                let mut cov_st = 0;
                let mut cov_en = 0;
                let mut cov = 0;

                for m in &b {
                    let st1 = std::cmp::max(st0, m.st);
                    let en1 = std::cmp::min(en0, m.en);
                    if st1 > cov_en {
                        cov += cov_en - cov_st;
                        cov_st = st1;
                        cov_en = en1;
                    } else {
                        cov_en = std::cmp::max(cov_en, en1);
                    }
                }
                cov += cov_en - cov_st;
                (st0, en0, b.len(), cov)
            }
            None => (st0, en0, 0, 0),
        };
        writeln!(stdout, "{}\t{}\t{}\t{}\t{}", chrom, start, end, len, cov)?;
        line.clear();
        b.clear();
        reader
            .read_line(&mut line)
            .expect("Failed to parse record.");
    }

    Ok(())
}
