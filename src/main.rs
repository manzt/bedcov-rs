use std::collections::HashMap;
use std::io::prelude::*;
use std::{env, fs, io};

#[derive(Debug, Clone, Copy)]
struct Interval {
    st: i32,
    en: i32,
    max: i32,
}

#[derive(Debug)]
struct IITree {
    a: Vec<Interval>,
    stack: [(usize, usize, usize); 64],
}

impl IITree {
    fn new() -> Self {
        IITree {
            a: Vec::new(),
            stack: [(0, 0, 0); 64],
        }
    }

    fn add(&mut self, st: i32, en: i32, max: i32) {
        self.a.push(Interval { st, en, max })
    }

    fn index(&mut self) -> usize {
        let ref mut a = self.a;
        a.sort_unstable_by_key(|k| k.st);

        let n = a.len();
        let mut last = 0;
        let mut last_i = 1;
        (0..n).step_by(2).for_each(|i| {
            last = a[i].en;
            a[i].max = a[i].en;
            last_i = i;
        });

        let mut k = 1;
        while 1 << k <= a.len() {
            let step = 1 << (k + 1);
            let i0 = (1 << k) - 1;

            for i in (i0..n).step_by(step) {
                let x = 1 << (k - 1);
                let max = std::cmp::max(a[i].en, a[i - x].max);
                let e = if i + x < a.len() { a[i + x].max } else { last };
                let max = std::cmp::max(e, max);
                a[i].max = max;
            }

            last_i = if (last_i >> k & 1) != 0 {
                last_i - (1 << (k - 1))
            } else {
                last_i + (1 << (k - 1))
            };

            if last_i < a.len() {
                last = std::cmp::max(last, a[last_i].max);
            }

            k += 1;
        }

        k - 1
    }

    fn overlap(&mut self, st: i32, en: i32, b: &mut Vec<Interval>) {
        let mut h = 0;
        while 1 << h <= self.a.len() {
            h += 1;
        }
        h -= 1;
        self.stack[0] = ((1 << h) - 1, h, 0);
        let mut t = 1;
        while t > 0 {
            t -= 1;
            let (x, h, w) = self.stack[t];
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
                self.stack[t] = (x, h, 1);
                t += 1;
                let y = x - (1 << (h - 1));
                if y >= self.a.len() || self.a[y].max > st {
                    self.stack[t] = (y, h - 1, 0);
                }
                t += 1;
            } else if x < self.a.len() && self.a[x].st < en {
                if st < self.a[x].en {
                    b.push(self.a[x]);
                }
                self.stack[t] = (x + (1 << (h - 1)), h - 1, 0);
                t += 1;
            }
        }
    }
}

fn parse_line(line: &str) -> Option<(&str, i32, i32)> {
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
    for line in content.lines() {
        let (chrom, start, end) = parse_line(line).expect("Failed to parse BED row.");
        let tree = map.entry(chrom).or_insert(IITree::new());
        tree.add(start, end, end);
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
