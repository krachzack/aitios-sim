#![feature(test)]

extern crate aitios_geom;
extern crate fixtures;
extern crate test;

use aitios_geom::{TangentSpace, Vec3};
use fixtures::venus::make_tracer;

#[bench]
fn trace_straight(b: &mut test::Bencher) {
    let tracer = make_tracer();
    let top = Vec3::new(0.0, 100.0, 0.0);
    let down = Vec3::new(0.0, -1.0, 0.0);

    b.iter(|| tracer.trace_straight(top, down).unwrap())
}

#[bench]
fn trace_parabolic(b: &mut test::Bencher) {
    let tracer = make_tracer();
    let top = Vec3::new(0.0, 100.0, 0.0);
    let up = Vec3::new(0.0, 1.0, 0.0);
    let down = -up;
    let hit = tracer.trace_straight(top, down).unwrap();
    let ontop = hit.intersection_point;
    let parabola_height = 1.0;

    b.iter(|| tracer.trace_parabolic(ontop, up, parabola_height).unwrap())
}

#[bench]
fn trace_flow(b: &mut test::Bencher) {
    let tracer = make_tracer();
    let top = Vec3::new(0.0, 100.0, 0.0);
    let up = Vec3::new(0.0, 1.0, 0.0);
    let down = -up;
    let hit = tracer.trace_straight(top, down).unwrap();
    let normal = hit.triangle.normal();
    let tangent = hit.triangle.tangent();
    let ontop = hit.intersection_point + 0.000001 * normal;

    b.iter(|| tracer.trace_flow(ontop, normal, tangent, 0.05).unwrap());
}
