#![feature(test)]

extern crate test;
extern crate fixtures;

use fixtures::venus::make_simulation;

/*#[bench]
fn simulate_iteration_with_1_ton(b: &mut test::Bencher) {
    let mut simulation = make_simulation(1);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}

#[bench]
fn simulate_iteration_with_1000_tons(b: &mut test::Bencher) {
    let mut simulation = make_simulation(1000);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}*/

#[bench]
fn simulate_iteration_with_10_000_tons(b: &mut test::Bencher) {
    let mut simulation = make_simulation(10000);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}

#[bench]
fn simulate_iteration_with_10_000_tons_parallelism(b: &mut test::Bencher) {
    let mut simulation = make_simulation(10000);

    b.iter(move || {
        simulation.run_fast();
        simulation.surfel_count()
    })
}

/*#[bench]
fn simulate_100_000(b: &mut test::Bencher) {
    let mut simulation = make_simulation(100000);

    b.iter(move || {
        simulation.run();
        simulation.surfel_count()
    })
}*/
