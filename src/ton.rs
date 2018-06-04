
use geom::{Vec3, Interpolation, TupleTriangle, Vertex, Position, TangentSpace};
use sampling::{TriangleBins, UnitSphere, UnitHemisphere, Uniform};
use scene::{Entity, Mesh};
use std::f32::EPSILON;
use std::iter;
use std::ops::Deref;

#[derive(Debug, Clone)]
pub struct Ton {
    /// Probability of moving further in a straight line
    pub p_straight: f32,
    /// Probability of moving further in a piecewise approximated
    /// parabolic trajectory
    pub p_parabolic: f32,
    /// Probability of moving tangently
    pub p_flow: f32,
    /// Determines the radius around a ton where it interacts with surface elements.
    pub interaction_radius: f32,
    /// Determines the height of a vertical bounce
    pub parabola_height: f32,
    /// Distance of a flow event
    pub flow_distance: f32,
    /// Flow direction calculation method
    pub flow_direction: FlowDirection,
    /// Amount of substances currently being carried by this ton
    pub substances: Vec<f32>,
    /// Factor by which the gammaton picks up material from surfels
    pub pickup_rates: Vec<f32>
}

/// Determines method for determination of flow direction based on a hit.
#[derive(Debug, Clone)]
pub enum FlowDirection {
    /// Flow direction is determined by projecting incoming direction of gammaton onto
    /// tangential plane.
    Incident,
    /// Projects the associated normalized vector onto the tangential plane to obtain
    /// flow direction.
    Static(Vec3)
}

// TODO the sampling should be stratified, e.g. by subdividing the possible directions into patches and ensuring every one gets its turn
enum Shape {
    /// A point source shooting equally in all directions
    Point { position: Vec3 },
    /// A hemispherical source aligned with the y axis shooting inward.
    Hemisphere {
        /// The center of the bottom disk of the hemisphere
        center: Vec3,
        /// Distance from the center for ray origins
        radius: f32
    },
    /// Shoots from the given mesh in interpolated normal direction
    Mesh { triangles: TriangleBins<TupleTriangle<Vertex>>, diffuse: bool }
}

pub struct TonSource {
    /// Emission shape
    shape: Shape,
    proto_ton: Ton,
    emission_count: usize
}

pub struct TonSourceBuilder {
    source: TonSource
}

pub struct TonEmission {
    pub origin: Vec3,
    pub direction: Vec3,
    pub ton: Ton
}

impl TonSource {

    fn emit_one(&self) -> TonEmission {
        let ton = self.proto_ton.clone();
        let (origin, direction) = match &self.shape {
            &Shape::Point { position } => (
                position.clone(),
                // Random position on the unit sphere
                UnitSphere.uniform()
            ),
            &Shape::Hemisphere { center, radius } => {
                let unit = UnitHemisphere::PosZ.uniform();
                let origin = center + radius * unit;
                // REVIEW wait, should they really all be flying towards the center?
                let direction = -unit;
                (origin, direction)
            },
            &Shape::Mesh { ref triangles, diffuse } => {
                // Interpolate a vertex on a random position on a randomly selected triangle (weighted by area)
                let tri = triangles.sample();
                let vtx = tri.interpolate_at(
                    tri.uniform(),
                    |v| v.clone()
                );

                let direction = if diffuse {
                    tri.tangent_to_world_matrix() * UnitHemisphere::PosZ.uniform()
                } else {
                    vtx.normal
                };

                let origin = vtx.position + direction * EPSILON;
                (origin, direction)
            }
        };

        TonEmission { origin, direction, ton }
    }

    pub fn emit<'a>(&'a self) -> impl Iterator<Item = TonEmission> + 'a {
        iter::repeat(self)
            .take(self.emission_count)
            .map(Self::emit_one)
    }

    pub fn emission_count(&self) -> usize {
        self.emission_count
    }
}

impl TonSourceBuilder {
    pub fn new() -> TonSourceBuilder {
        TonSourceBuilder {
            source: TonSource {
                emission_count: 10000,
                shape: Shape::Point { position: Vec3::new(0.0, 0.0, 0.0) },
                proto_ton: Ton {
                    p_straight: 0.0,
                    p_parabolic: 0.0,
                    p_flow: 0.0,
                    substances: Vec::new(),
                    interaction_radius: 0.1,
                    parabola_height: 0.05,
                    flow_distance: 0.02,
                    flow_direction: FlowDirection::Incident,
                    pickup_rates: Vec::new()
                }
            }
        }
    }

    pub fn point_shaped(mut self, pos_x: f32, pos_y: f32, pos_z: f32) -> TonSourceBuilder {
        self.source.shape = Shape::Point { position: Vec3::new(pos_x, pos_y, pos_z) };
        self
    }

    pub fn hemisphere_shaped(mut self, center: Vec3, radius: f32) -> TonSourceBuilder {
        self.source.shape = Shape::Hemisphere { center, radius };
        self
    }

    pub fn entity_shaped(self, entity: &Entity, diffuse: bool) -> TonSourceBuilder {
        self.mesh_shaped(&entity.mesh, diffuse)
    }

    pub fn mesh_shaped<'a, T, M, V>(mut self, mesh: &'a T, diffuse: bool) -> TonSourceBuilder
        where T : Deref<Target = M>,
            M : Mesh<'a, Vertex = V> + 'a,
            V : Position,
            TriangleBins<TupleTriangle<Vertex>>: iter::FromIterator<TupleTriangle<V>>
    {
        self.source.shape = Shape::Mesh {
            triangles: mesh.triangles()
                .collect(),
            diffuse
        };

        self
    }

    pub fn emission_count(mut self, emission_count: usize) -> TonSourceBuilder {
        self.source.emission_count = emission_count;
        self
    }

    pub fn p_straight(mut self, p_straight: f32) -> TonSourceBuilder {
        self.source.proto_ton.p_straight = p_straight;
        self
    }

    pub fn p_parabolic(mut self, p_parabolic: f32) -> TonSourceBuilder {
        self.source.proto_ton.p_parabolic = p_parabolic;
        self
    }

    pub fn p_flow(mut self, p_flow: f32) -> TonSourceBuilder {
        self.source.proto_ton.p_flow = p_flow;
        self
    }

    pub fn substances(mut self, substances: &Vec<f32>) -> TonSourceBuilder {
        self.source.proto_ton.substances = substances.clone();
        self
    }

    pub fn pickup_rates<R : IntoIterator<Item = f32>> (mut self, pickup_rates: R) -> TonSourceBuilder {
        self.source.proto_ton.pickup_rates = pickup_rates.into_iter().collect();
        self
    }

    pub fn interaction_radius(mut self, interaction_radius: f32) -> TonSourceBuilder {
        self.source.proto_ton.interaction_radius = interaction_radius;
        self
    }

    pub fn parabola_height(mut self, parabola_height: f32) -> TonSourceBuilder {
        self.source.proto_ton.parabola_height = parabola_height;
        self
    }

    pub fn flow_distance(mut self, flow_distance: f32) -> TonSourceBuilder {
        self.source.proto_ton.flow_distance = flow_distance;
        self
    }

    pub fn flow_direction_static(mut self, flow_direction: Vec3) -> TonSourceBuilder {
        self.source.proto_ton.flow_direction = FlowDirection::Static(flow_direction);
        self
    }

    pub fn build(self) -> TonSource {
        assert_eq!(
            self.source.proto_ton.pickup_rates.len(),
            self.source.proto_ton.substances.len(),
            "Pickup rates and initial substance concentrations have unequal lengths"
        );

        self.source
    }
}

#[cfg(test)]
mod test {
    extern crate aitios_asset;

    use super::*;

    #[test]
    fn test_shoot_from_mesh() {
        let entities = aitios_asset::obj::load("test-scenes/buddha-scene-ton-source-mesh/buddha-scene-ton-source-sun.obj")
            .unwrap();

        assert_eq!(1, entities.len());

        let src = TonSourceBuilder::new()
            .p_flow(0.2)
            .emission_count(10)
            .entity_shaped(&entities[0], false)
            .build();

        assert_eq!(src.emit().count(), 10);
        assert!(src.emit().all(|TonEmission { ton, origin, direction }| ton.p_flow == 0.2 && origin.y > 0.1 && direction.y < 0.0));

        let src = TonSourceBuilder::new()
            .p_flow(0.2)
            .emission_count(10)
            .mesh_shaped(&entities[0].mesh, false)
            .build();

        assert_eq!(src.emit().count(), 10);
        assert!(src.emit().all(|TonEmission { ton, origin, direction }| ton.p_flow == 0.2 && origin.y > 0.1 && direction.y < 0.0));
    }
}
