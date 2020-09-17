use linereader::LineReader;
use std::collections::HashMap;
use std::io::prelude::*;
use std::{env, fs, io};

#[derive(Debug)]
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

    fn index(&mut self) {
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
    }

    fn overlap(&mut self, st: i32, en: i32) -> (i32, i32, i32, i32) {
        let mut cov_st = 0;
        let mut cov_en = 0;
        let mut cov = 0;
        let mut b = 0;

        let mut calc = |i: &Interval| {
            let st1 = std::cmp::max(st, i.st);
            let en1 = std::cmp::min(en, i.en);
            if st1 > cov_en {
                cov += cov_en - cov_st;
                cov_st = st1;
                cov_en = en1;
            } else {
                cov_en = std::cmp::max(cov_en, en1);
            }
            b += 1;
        };

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
                        calc(&self.a[i])
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
                    calc(&self.a[x])
                }
                self.stack[t] = (x + (1 << (h - 1)), h - 1, 0);
                t += 1;
            }
        }
        cov += cov_en - cov_st;
        (st, en, b, cov)
    }
}

#[inline]
fn parse_line(line: &[u8]) -> Option<(&str, i32, i32)> {
    let line = std::str::from_utf8(line).ok()?;
    let mut iter = line.trim().split("\t").take(3);
    let chrom = iter.next()?;
    let start = iter.next()?.parse().ok()?;
    let end = iter.next()?.parse().ok()?;
    Some((chrom, start, end))
}

fn main() -> io::Result<()> {
    let path = env::args().nth(1).expect("Missing input file.");
    let f = fs::File::open(path)?;
    let mut reader = LineReader::new(f);

    let mut map: HashMap<String, IITree> = HashMap::with_capacity(24);
    while let Some(line) = reader.next_line() {
        let line = line.expect("Read error.");
        let (chrom, start, end) = parse_line(&line).expect("Failed to parse BED row.");
        match map.get_mut(chrom) {
            Some(tree) => tree.add(start, end, end),
            None => {
                let mut tree = IITree::new();
                tree.add(start, end, end);
                map.insert(chrom.to_string(), tree);
            }
        }
    }

    map.values_mut().for_each(|tree| tree.index());

    let path = env::args().nth(2).expect("Missing file 2.");
    let f = fs::File::open(path)?;
    let mut reader = LineReader::new(f);

    let stdout = io::stdout();
    let mut stdout = stdout.lock();

    while let Some(line) = reader.next_line() {
        let line = line.expect("Read error.");
        let (chrom, st0, en0) = parse_line(&line).expect("failed to parse row.");
        let (start, end, len, cov) = match map.get_mut(chrom) {
            Some(tree) => tree.overlap(st0, en0),
            None => (st0, en0, 0, 0),
        };
        writeln!(stdout, "{}\t{}\t{}\t{}\t{}", chrom, start, end, len, cov)?;
    }

    Ok(())
}
