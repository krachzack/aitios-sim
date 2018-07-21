#![feature(test)]

extern crate fixtures;
extern crate test;

use fixtures::venus::make_simulation;

#[bench]
fn simulate_iteration_with_1_000_tons(b: &mut test::Bencher) {
    let mut simulation = make_simulation(1000);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}

#[bench]
fn simulate_iteration_with_10_000_tons(b: &mut test::Bencher) {
    let mut simulation = make_simulation(10000);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}

#[ignore]
#[bench]
fn simulate_iteration_with_100_000_tons(b: &mut test::Bencher) {
    let mut simulation = make_simulation(100000);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}

#[ignore]
#[bench]
fn simulate_iteration_with_1_000_000_tons(b: &mut test::Bencher) {
    let mut simulation = make_simulation(100000);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}
