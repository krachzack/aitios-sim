
use geom::{Vec3, Interpolation, TupleTriangle, Vertex, Position, TangentSpace};
use sampling::{TriangleBins, UnitSphere, UnitHemisphere, Uniform};
use scene::{Entity, Mesh};
use std::f32::EPSILON;
use std::iter;

#[derive(Debug)]
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
    /// Amount of substances currently being carried by this ton
    pub substances: Vec<f32>,
    /// Factor by which the gammaton picks up material from surfels
    pub pickup_rates: Vec<f32>
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
    /// Probability of moving further in a straight line for tons emitted by this source
    p_straight: f32,
    /// Probability of moving further in a piecewise approximated for tons emitted by this source
    /// parabolic trajectory
    p_parabolic: f32,
    /// Probability of moving tangently for tons emitted by this source
    p_flow: f32,
    /// Determines the radius around a ton where it interacts with surface elements.
    interaction_radius: f32,
    /// Determines the height of a vertical bounce
    parabola_height: f32,
    /// Distance of a flow event
    flow_distance: f32,
    /// Amount of substances initially carried by tons emitted by this source
    substances: Vec<f32>,
    emission_count: usize,
    pickup_rates: Vec<f32>
}

pub struct TonSourceBuilder {
    /// Emission shape
    shape: Shape,
    /// Probability of moving further in a straight line for tons emitted by this source
    p_straight: f32,
    /// Probability of moving further in a piecewise approximated for tons emitted by this source
    /// parabolic trajectory
    p_parabolic: f32,
    /// Probability of moving tangently for tons emitted by this source
    p_flow: f32,
    /// Amount of substances initially carried by tons emitted by this source
    substances: Vec<f32>,
    emission_count: usize,
    pickup_rates: Vec<f32>,
    /// Determines the radius around a ton where it interacts with surface elements.
    interaction_radius: f32,
    /// Determines the height of a vertical bounce
    parabola_height: f32,
    /// Distance of a flow event
    flow_distance: f32,
}

pub struct TonEmission {
    pub origin: Vec3,
    pub direction: Vec3,
    pub ton: Ton
}

impl TonSource {

    pub fn emit_one(&self) -> TonEmission {
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

        let ton = Ton {
            p_straight: self.p_straight,
            p_parabolic: self.p_parabolic,
            p_flow: self.p_flow,
            interaction_radius: self.interaction_radius,
            parabola_height: self.parabola_height,
            flow_distance: self.flow_distance,
            substances: self.substances.clone(),
            pickup_rates: self.pickup_rates.clone()
        };

        TonEmission { origin, direction, ton }
    }

    pub fn emit<'a>(&'a self) -> Box<Iterator<Item=TonEmission> + 'a> {
        Box::new(iter::repeat(self)
            .map(Self::emit_one)
            .take(self.emission_count as usize))
    }

    pub fn emission_count(&self) -> usize {
        self.emission_count
    }
}

impl TonSourceBuilder {
    pub fn new() -> TonSourceBuilder {
        TonSourceBuilder {
            p_straight: 0.0,
            p_parabolic: 0.0,
            p_flow: 0.0,
            substances: Vec::new(),
            shape: Shape::Point { position: Vec3::new(0.0, 0.0, 0.0) },
            emission_count: 10000,
            interaction_radius: 0.1,
            parabola_height: 0.05,
            flow_distance: 0.02,
            pickup_rates: Vec::new()
        }
    }

    pub fn p_straight(mut self, p_straight: f32) -> TonSourceBuilder {
        self.p_straight = p_straight;
        self
    }

    #[allow(dead_code)]
    pub fn p_parabolic(mut self, p_parabolic: f32) -> TonSourceBuilder {
        self.p_parabolic = p_parabolic;
        self
    }

    #[allow(dead_code)]
    pub fn p_flow(mut self, p_flow: f32) -> TonSourceBuilder {
        self.p_flow = p_flow;
        self
    }

    pub fn substances(mut self, substances: &Vec<f32>) -> TonSourceBuilder {
        self.substances = substances.clone();
        self
    }

    pub fn point_shaped(mut self, pos_x: f32, pos_y: f32, pos_z: f32) -> TonSourceBuilder {
        self.shape = Shape::Point { position: Vec3::new(pos_x, pos_y, pos_z) };
        self
    }

    pub fn hemisphere_shaped(mut self, center: Vec3, radius: f32) -> TonSourceBuilder {
        self.shape = Shape::Hemisphere { center, radius };
        self
    }

    pub fn entity_shaped(self, entity: &Entity, diffuse: bool) -> TonSourceBuilder {
        self.mesh_shaped(&entity.mesh, diffuse)
    }

    pub fn mesh_shaped<'a, M, V>(mut self, mesh: &'a M, diffuse: bool) -> TonSourceBuilder
        where M : Mesh<'a, Vertex = V>,
            V : Position,
            TriangleBins<TupleTriangle<Vertex>>: iter::FromIterator<TupleTriangle<V>>
    {
        self.shape = Shape::Mesh {
            triangles: mesh.triangles()
                .collect(),
            diffuse
        };

        self
    }

    pub fn emission_count(mut self, emission_count: usize) -> TonSourceBuilder {
        self.emission_count = emission_count;
        self
    }

    pub fn interaction_radius(mut self, interaction_radius: f32) -> TonSourceBuilder {
        self.interaction_radius = interaction_radius;
        self
    }

    pub fn parabola_height(mut self, parabola_height: f32) -> TonSourceBuilder {
        self.parabola_height = parabola_height;
        self
    }

    pub fn flow_distance(mut self, flow_distance: f32) -> TonSourceBuilder {
        self.flow_distance = flow_distance;
        self
    }

    pub fn pickup_rates<R : IntoIterator<Item = f32>> (mut self, pickup_rates: R) -> TonSourceBuilder {
        self.pickup_rates = pickup_rates.into_iter().collect();
        self
    }

    pub fn build(self) -> TonSource {
        assert_eq!(self.pickup_rates.len(), self.substances.len());

        TonSource {
            shape: self.shape,
            p_straight: self.p_straight,
            p_parabolic: self.p_parabolic,
            p_flow: self.p_flow,
            interaction_radius: self.interaction_radius,
            parabola_height: self.parabola_height,
            flow_distance: self.flow_distance,
            substances: self.substances,
            emission_count: self.emission_count,
            pickup_rates: self.pickup_rates
        }
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
