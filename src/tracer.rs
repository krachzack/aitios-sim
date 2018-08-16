use geom::prelude::*;
use geom::{TupleTriangle, Vec3, Vertex};
use spatial::Octree;
#[cfg(feature = "debug_tracing")]
use std::cell::RefCell;
use std::f32::INFINITY;

#[cfg(feature = "debug_tracing")]
enum TracingEvent {
    Straight(Vec3, Vec3),
    Parabolic(Vec3, Vec3),
    Flow(Vec3, Vec3),
}
#[cfg(feature = "debug_tracing")]
const MAX_TRACING_EVENT_COUNT: usize = 600;
/// Delta to move before checking for intersections. This avoids rays
/// intersecting the triangle they originated from due to floating point
/// imprecision.
const SELF_INTERSECTION_EPSILON: f32 = 0.0001;
const FLOW_ADHESIVENESS: f32 = 1.1; // Allow downward motion of flow to be 10% longer than to be expected on a flat surface, so it can flow upward a little bit

pub struct Tracer {
    geometry: Octree<TupleTriangle<Vertex>>,
    gravity_direction: Vec3,
    #[cfg(feature = "debug_tracing")]
    first_tracing_events: RefCell<Vec<TracingEvent>>,
}

#[derive(Debug)]
pub struct Hit<'a> {
    pub intersection_point: Vec3,
    pub incoming_direction: Vec3,
    pub triangle: &'a TupleTriangle<Vertex>,
}

impl Tracer {
    pub fn new<I>(triangles: I) -> Self
    where
        I: IntoIterator<Item = TupleTriangle<Vertex>>,
    {
        Tracer {
            geometry: triangles.into_iter().collect(),
            gravity_direction: Vec3::new(0.0, -1.0, 0.0),
            #[cfg(feature = "debug_tracing")]
            first_tracing_events: RefCell::new(Vec::new()),
        }
    }

    pub fn trace_straight(&self, from: Vec3, direction: Vec3) -> Option<Hit> {
        let from = from + direction * 0.0000001; // avoid self-intersection
        self.geometry
            .ray_intersection_target_and_parameter(from, direction)
            .map(|(hit_tri, t)| {
                let intersection_point = from + t * direction;

                #[cfg(feature = "debug_tracing")]
                self.debug_straight(from, intersection_point);

                Hit {
                    intersection_point,
                    incoming_direction: direction,
                    triangle: hit_tri,
                }
            })
    }

    /// `upward_parabola_height` indicates maximum height of a bounce in the special case it is flung
    /// straight up and gravity pointing straight down. Most parabolic paths will not be straight up and
    /// have less height, might even be pointing downward from the beginning.
    pub fn trace_parabolic(
        &self,
        from: Vec3,
        direction: Vec3,
        upward_parabola_height: f32,
    ) -> Option<Hit> {
        let gravity_mag = 9.81_f32;
        let timestep = 1.0 / 60.0; // 0.0333333 seconds, more is more exact but slower

        let gravity_acceleration = gravity_mag * self.gravity_direction;
        let takeoff_velocity_mag = (2.0 * gravity_mag * upward_parabola_height).sqrt();
        let mut velocity = takeoff_velocity_mag * direction;
        let mut position = from + direction * SELF_INTERSECTION_EPSILON;

        let scene_bounds = {
            let mut bounds = self.geometry.bounds();
            bounds.max.y = INFINITY; // ignore if out of bounds in positive y direction since gravity will eventually pull it downward
            bounds
        };

        // REVIEW This should take into account that the source could be outside the scene bounds.
        //        At emission time, the tracing will be always straight though, so it should not come to that
        while scene_bounds.is_point_inside(position) {
            velocity += gravity_acceleration * timestep;

            let spatial_delta = velocity * timestep;
            let dist = spatial_delta.magnitude();
            let direction = spatial_delta / dist;

            if let Some((hit_tri, t)) = self
                .geometry
                .line_segment_intersection_target_and_parameter(position, direction, dist)
            {
                let intersection_point = position + t * direction;

                #[cfg(feature = "debug_tracing")]
                self.debug_parabolic(position, intersection_point);

                return Some(Hit {
                    intersection_point,
                    incoming_direction: direction,
                    triangle: hit_tri,
                });
            } else {
                // No intersection, safe to move particle without penetrating objects
                position += spatial_delta;

                #[cfg(feature = "debug_tracing")]
                self.debug_parabolic(position - spatial_delta, position);
            }
        }

        None
    }

    /// From should be on a triangle, direction should be aligned with tangential plane of triangle
    /// Up is in normal direction.
    /// flow_distance is offset in tangential direction before interacting again.
    pub fn trace_flow(
        &self,
        from: Vec3,
        up: Vec3,
        tangential_direction: Vec3,
        flow_distance: f32,
    ) -> Option<Hit> {
        // Upward epsilon is chosen with a fixed angle
        let upward_epsilon = flow_distance * 1.3;

        // FIXME
        // in the thesis it sounded better to first do
        // tangential, and the second segment either straight
        // down or gravity

        // First a little bias to avoid self-intersection
        let from = from + SELF_INTERSECTION_EPSILON * up;

        // Expected flow target on a hypothetical flat surface
        let to = from + tangential_direction * flow_distance;

        // First, move up for flow
        // This should normally not intersect anything, but
        // provides an origin for the downward raycast.
        // If it does intersect something, we are in some sort
        // of cavity. Count as flow target even though not tangential.
        let upward_hit = self
            .geometry
            .line_segment_intersection_target_and_parameter(from, up, upward_epsilon);

        if let Some((hit_tri, t)) = upward_hit {
            let intersection_point = from + t * up;

            #[cfg(feature = "debug_tracing")]
            self.debug_flow(from, intersection_point);

            return Some(Hit {
                intersection_point,
                incoming_direction: up,
                triangle: hit_tri,
            });
        }

        // Now move from above to projected intersection location
        let atop = from + upward_epsilon * up;
        let dir = (to - atop).normalize();
        // On a flat surface, the ray should intersect
        // again in sqrt(expected_dist_sqr) distance from
        // the new from offset to the top
        // TODO cache this, sqrt is expensive
        let expected_dist_sqr = upward_epsilon * upward_epsilon + flow_distance * flow_distance;
        let expected_dist = expected_dist_sqr.sqrt();

        #[cfg(feature = "debug_tracing")]
        self.debug_flow(from, atop);
        let tangential_hit = self
            .geometry
            .line_segment_intersection_target_and_parameter(
                atop,
                dir,
                expected_dist + FLOW_ADHESIVENESS,
            );

        if let Some((hit_tri, t)) = tangential_hit {
            let intersection_point = atop + t * dir;

            #[cfg(feature = "debug_tracing")]
            self.debug_flow(atop, intersection_point);

            return Some(Hit {
                intersection_point,
                incoming_direction: dir,
                triangle: hit_tri,
            });
        }

        #[cfg(feature = "debug_tracing")]
        self.debug_flow(atop, to);

        // Not back on the surface yet, surface must have
        // convex local neighbourhood. Follow gravity.
        let from = atop + dir * (expected_dist + FLOW_ADHESIVENESS);
        let secondary_target = self
            .geometry
            .ray_intersection_target_and_parameter(from, self.gravity_direction);

        if let Some((hit_tri, t)) = secondary_target {
            let intersection_point = from + self.gravity_direction * t;

            #[cfg(feature = "debug_tracing")]
            self.debug_flow(from, intersection_point);

            return Some(Hit {
                intersection_point,
                incoming_direction: self.gravity_direction,
                triangle: hit_tri,
            });
        }

        None
    }

    #[cfg(feature = "debug_tracing")]
    fn debug_push(&self, evt: TracingEvent) {
        let len = self.first_tracing_events.borrow().len();

        if len == MAX_TRACING_EVENT_COUNT {
            return;
        } else if len == MAX_TRACING_EVENT_COUNT - 1 {
            {
                let mut evts = self.first_tracing_events.borrow_mut();
                evts.push(evt);
            }
            self.dump_debug_tracing_events();
        } else {
            let mut evts = self.first_tracing_events.borrow_mut();
            evts.push(evt);
        }
    }

    #[cfg(feature = "debug_tracing")]
    fn debug_straight(&self, from: Vec3, to: Vec3) {
        self.debug_push(TracingEvent::Straight(from, to));
    }

    #[cfg(feature = "debug_tracing")]
    fn debug_parabolic(&self, from: Vec3, to: Vec3) {
        self.debug_push(TracingEvent::Parabolic(from, to));
    }

    #[cfg(feature = "debug_tracing")]
    fn debug_flow(&self, from: Vec3, to: Vec3) {
        self.debug_push(TracingEvent::Flow(from, to));
    }

    #[cfg(feature = "debug_tracing")]
    fn dump_debug_tracing_events(&self) {
        use std::fs::File;
        use std::io::Write;
        use std::path::PathBuf;

        let obj_path = PathBuf::from("debug_tracing.obj");
        /*if obj_path.is_file() {
            return; // Already dumped
        }*/

        let mtl_path = PathBuf::from("debug_tracing.mtl");
        if !mtl_path.is_file() {
            let mut mtl = File::create(&mtl_path).unwrap();

            // Straight has green diffuse color
            mtl.write("\nnewmtl straight\n".as_bytes()).unwrap();
            mtl.write("Ka 1.0 1.0 1.0\n".as_bytes()).unwrap();
            mtl.write("Kd 0.0 1.0 0.0\n".as_bytes()).unwrap();
            mtl.write("Ks 1.0 1.0 1.0\n".as_bytes()).unwrap();
            mtl.write("illum 2\n".as_bytes()).unwrap();

            // Parabolic has blue diffuse color
            mtl.write("\nnewmtl parabolic\n".as_bytes()).unwrap();
            mtl.write("Ka 1.0 1.0 1.0\n".as_bytes()).unwrap();
            mtl.write("Kd 0.0 0.0 1.0\n".as_bytes()).unwrap();
            mtl.write("Ks 1.0 1.0 1.0\n".as_bytes()).unwrap();
            mtl.write("illum 2\n".as_bytes()).unwrap();

            // Parabolic has red diffuse color
            mtl.write("\nnewmtl flow\n".as_bytes()).unwrap();
            mtl.write("Ka 1.0 1.0 1.0\n".as_bytes()).unwrap();
            mtl.write("Kd 1.0 0.0 0.0\n".as_bytes()).unwrap();
            mtl.write("Ks 1.0 1.0 1.0\n".as_bytes()).unwrap();
            mtl.write("illum 2\n".as_bytes()).unwrap();
        }

        let mut obj = File::create(&obj_path).unwrap();
        obj.write("mtllib debug_tracing.mtl\n".as_bytes()).unwrap();
        let mut idx = 1;

        let evts = self.first_tracing_events.borrow();
        let mut last_mtl = "";
        for evt in evts.iter() {
            let (mat, from, to) = match evt {
                &TracingEvent::Straight(from, to) => ("straight", from, to),
                &TracingEvent::Parabolic(from, to) => ("parabolic", from, to),
                &TracingEvent::Flow(from, to) => ("flow", from, to),
            };

            if mat != last_mtl {
                obj.write(format!("o {}\n", mat).as_bytes()).unwrap();
                obj.write(format!("usemtl {}\n", mat).as_bytes()).unwrap();
                last_mtl = mat;
            }

            obj.write(format!("v {} {} {}\n", from.x, from.y, from.z).as_bytes())
                .unwrap();
            obj.write(format!("v {} {} {}\n", to.x, to.y, to.z).as_bytes())
                .unwrap();
            obj.write(format!("l {} {}\n", idx, idx + 1).as_bytes())
                .unwrap();

            idx += 2;
        }
    }
}

#[cfg(test)]
mod test {
    extern crate aitios_asset;

    use super::*;
    use geom::{Position, TangentSpace, TupleTriangle as Tri, Vec2, Vec3, Vertex};
    use scene::Mesh;

    #[test]
    fn test_straight_tracing() {
        // This is the top part of an icosphere, will try to hit it by shooting from the origin
        let entities = aitios_asset::obj::load(
            "test-scenes/buddha-scene-ton-source-mesh/buddha-scene-ton-source-sun.obj",
        ).unwrap();

        // Will try to hit first vertex in first entity
        let known_vertex = entities[0].mesh.vertices().next().unwrap().position();

        let tracer = Tracer::new(entities.iter().flat_map(|ent| ent.mesh.triangles()));

        let origin = Vec3::new(0.0, 0.1, 0.2);
        let direction = known_vertex - origin;

        // When shooting from the origin, should hit it
        let hit = tracer.trace_straight(origin, direction);
        assert!(hit.is_some(), "Expected to hit known vertex");
        assert_ulps_eq!(hit.unwrap().intersection_point, known_vertex);

        // In the other direction, it should be a miss
        let miss = tracer.trace_straight(origin, -direction);
        assert!(
            miss.is_none(),
            "Expected to hit nothing when shooting in inverted direction"
        );
    }

    #[test]
    fn test_straight_up_parabola() {
        // This is the top part of an icosphere, will try to hit it by shooting a vertex
        // from exactly below and a given parabola height
        let entities = aitios_asset::obj::load(
            "test-scenes/buddha-scene-ton-source-mesh/buddha-scene-ton-source-sun.obj",
        ).unwrap();

        let known_vertex = entities[0].mesh.vertices().next().unwrap().position();

        let tracer = Tracer::new(entities.iter().flat_map(|ent| ent.mesh.triangles()));

        let parabola_height = 10.0;
        let bias = 0.2; // The bias is necessary because euler integration is quite inexact, especially with f32

        let origin = known_vertex + Vec3::new(0.0, -parabola_height + bias, 0.0);

        let hit = tracer.trace_parabolic(origin, Vec3::new(0.0, 1.0, 0.0), parabola_height);
        assert!(hit.is_some(), "Expected to hit known vertex");
        assert_ulps_eq!(hit.unwrap().intersection_point, known_vertex);
    }

    #[test]
    fn test_flow_on_flat() {
        // 1x1 quad on the X/Z plane
        let tracer = Tracer::new(x_z_quad().into_iter());

        let origin = Vec3::new(0.0, 0.0, 0.0);
        let up = Vec3::new(0.0, 1.0, 0.0);
        let flow_direction = Vec3::new(0.0, 0.0, 1.0);

        let hit = tracer.trace_flow(origin, up, flow_direction, 0.4);
        assert!(hit.is_some());
        let hit = hit.unwrap();
        assert_relative_eq!(
            Vec3::new(0.0, 0.0, 0.4),
            hit.intersection_point,
            epsilon = 0.0001
        );

        let hit = tracer.trace_flow(origin, up, flow_direction, 0.8);
        assert!(hit.is_some());
        let hit = hit.unwrap();
        assert_relative_eq!(
            Vec3::new(0.0, 0.0, 0.8),
            hit.intersection_point,
            epsilon = 0.0001
        );

        let hit = tracer.trace_flow(origin, up, flow_direction, 0.9);
        assert!(hit.is_some());
        let hit = hit.unwrap();
        assert_relative_eq!(
            Vec3::new(0.0, 0.0, 0.9),
            hit.intersection_point,
            epsilon = 0.0001
        );

        let hit = tracer.trace_flow(origin, up, flow_direction, 0.99);
        assert!(hit.is_some());
        let hit = hit.unwrap();
        assert_relative_eq!(
            Vec3::new(0.0, 0.0, 0.99),
            hit.intersection_point,
            epsilon = 0.0001
        );

        let hit = tracer
            .trace_flow(origin, up, flow_direction, 1.5)
            .map(|h| h.intersection_point);
        assert_eq!(None, hit);

        let hit = tracer.trace_flow(origin, up, flow_direction, 1.05);
        assert!(hit.is_none());

        let hit = tracer.trace_flow(origin, up, flow_direction, 1.01);
        assert!(hit.is_none());
    }

    #[test]
    #[ignore]
    fn test_flow_on_concave() {
        unimplemented!();
    }

    #[test]
    #[ignore]
    fn test_flow_on_convex() {
        unimplemented!();
    }

    #[test]
    fn test_flow_on_flat_surface_with_offset() {
        // This is the top part of an icosphere, will try to hit it by shooting a vertex
        // from exactly below and a given parabola height
        let entities = aitios_asset::obj::load(
            "test-scenes/buddha-scene-ton-source-mesh/buddha-scene-ton-source-sun.obj",
        ).unwrap();

        let tracer = Tracer::new(entities.iter().flat_map(|ent| ent.mesh.triangles()));

        let some_tri = entities[0].mesh.triangles().next().unwrap();

        let centroid = some_tri.centroid();
        let normal = some_tri.normal();
        let (vertex_a, _, _) = some_tri.positions();

        let to_vertex_a_dist = vertex_a.distance(centroid);
        let flow_direction = (centroid - vertex_a) / to_vertex_a_dist;
        let expected_intersection_point = centroid;

        let hit = tracer.trace_flow(vertex_a, normal, flow_direction, to_vertex_a_dist);

        assert!(hit.is_some(), "Expected to hit known vertex");

        let hit = hit.unwrap().intersection_point;
        assert!(
            hit.distance(expected_intersection_point) < 0.0001,
            "Flow tracing ended up in unexpected place: {:?}\n expecting: {:?}",
            hit,
            expected_intersection_point
        );
    }

    fn x_z_quad() -> Vec<Tri<Vertex>> {
        let left_front = Vertex {
            position: Vec3::new(-1.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            texcoords: Vec2::new(0.0, 0.0),
        };

        let right_front = Vertex {
            position: Vec3::new(1.0, 0.0, 1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            texcoords: Vec2::new(1.0, 0.0),
        };

        let left_rear = Vertex {
            position: Vec3::new(-1.0, 0.0, -1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            texcoords: Vec2::new(0.0, 1.0),
        };

        let right_rear = Vertex {
            position: Vec3::new(1.0, 0.0, -1.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            texcoords: Vec2::new(1.0, 1.0),
        };

        vec![
            Tri(left_front.clone(), right_front, right_rear.clone()),
            Tri(left_front, right_rear, left_rear),
        ]
    }
}
