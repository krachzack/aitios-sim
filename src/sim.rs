use ::tracer::{Tracer, Hit};
use ::surfel_data::SurfelData;
use ::surfel_rule::SurfelRule;
use ::ton::{Ton, TonSource, TonEmission, FlowDirection};
use geom::{Vertex, TupleTriangle, TangentSpace};
use geom::prelude::*;
use surf;
use surf::Surfel;
use sampling::{Uniform, UnitHemisphere};
use rand;
use rand::Rng;

type Surface = surf::Surface<Surfel<Vertex, SurfelData>>;

pub struct Simulation {
    sources: Vec<TonSource>,
    tracer: Tracer,
    surface: Surface,
    /// Global surfel rules for all surfels
    surfel_rules: Vec<SurfelRule>
}

#[derive(PartialEq, Copy, Clone)]
enum MotionType {
    Straight,
    Parabolic,
    Flow,
    Settled
}

impl Simulation {
    pub fn new<I>(sources: Vec<TonSource>, triangles: I, surface: Surface, surfel_rules: Vec<SurfelRule>) -> Self
        where I : IntoIterator<Item = TupleTriangle<Vertex>>
    {
        Simulation {
            sources,
            surface,
            tracer: Tracer::new(triangles),
            surfel_rules
        }
    }

    // Advances the simulation by one iteration
    pub fn run(&mut self) {
        let sources = &self.sources;
        let tracer = &self.tracer;
        let rules = &self.surfel_rules;
        let surface = &mut self.surface;

        // TODO change reflectance properties based on gammaton map from last iteration
        //      to increase rusting rate in rusty areas

        sources.iter()
            .flat_map(TonSource::emit)
            .for_each(|e| Self::emit(surface, tracer, e));

        Self::perform_rules(surface, rules);
    }

    pub fn surface(&self) -> &Surface {
        &self.surface
    }

    pub fn surfel_count(&self) -> usize {
        self.surface.samples.len()
    }

    /// Amount of gammatons emitted from all sources each iteration
    pub fn emission_count(&self) -> usize {
        self.sources.iter()
            .map(|s| s.emission_count())
            .sum()
    }

    fn perform_rules(surf: &mut Surface, global_rules: &Vec<SurfelRule>) {
        // First the global rules
        for rule in global_rules {
            surf.samples.iter_mut()
                .for_each(|s| Self::perform_rule(&mut s.data_mut().substances, rule))
        }

        // Then the local ones
        surf.samples.iter_mut()
            .for_each(|s| {
                let s = s.data_mut();
                let substances = &mut s.substances;
                for rule in &s.rules {
                    Self::perform_rule(substances, rule);
                }
            });
    }

    fn perform_rule(substances: &mut Vec<f32>, rule: &SurfelRule) {
        // REVIEW should the substances be clamped?
        match rule {
            &SurfelRule::Deteriorate { substance_idx, factor } =>
                substances[substance_idx] = ((1.0 + factor) * substances[substance_idx]).max(0.0),

            &SurfelRule::Transfer { source_substance_idx, target_substance_idx, factor } => {
                let transport_amount = factor * substances[source_substance_idx];
                substances[source_substance_idx] = (substances[source_substance_idx] - transport_amount).max(0.0);
                substances[target_substance_idx] = (substances[target_substance_idx] + transport_amount).max(0.0);
            }
        }
    }

    fn emit(surf: &mut Surface, tracer: &Tracer, mut emission: TonEmission) {
        if let Some(hit) = tracer.trace_straight(emission.origin, emission.direction) {
            Self::interact(surf, tracer, &mut emission.ton, hit);
        }
    }

    fn interact(surf: &mut Surface, tracer: &Tracer, ton: &mut Ton, hit: Hit) {
        let mut surfel_idxs = surf.find_within_sphere_indexes(
            hit.intersection_point,
            ton.interaction_radius
        );

        // Throw out all surfels where normals are rotated by more than 90°
        // relative to the normal of the hit triangle.
        // This aims to minimize surfels from the other side of thin surfaces,
        // being affected from hits to the other side.
        // Depending on the interaction radius and the complexity of the
        // surface in the interaction radius range, bleeding may still occur.
        let hit_tri_normal = hit.triangle.normal();
        surfel_idxs.retain(|&i| {
            let surfel_normal  = surf.samples[i].vertex().normal;
            hit_tri_normal.dot(surfel_normal) > 0.0
        });

        if surfel_idxs.len() == 0 {
            warn!("Ton hit a surface but did not interact with any surfels, try higher interaction radius, interacting with nearest surfel instead.");
            surfel_idxs.push(surf.nearest_idx(hit.intersection_point));
        }

        let next_motion_type = Self::select_motion_type(ton);
        let count_weight = (surfel_idxs.len() as f32).recip();

        if next_motion_type == MotionType::Settled {
            // settle
            for idx in surfel_idxs {
                Self::deposit(ton, surf.samples[idx].data_mut(), count_weight);
            }
        } else {
            // bounce
            Self::deteriorate(ton, surf.samples[surfel_idxs[0]].data());

            for idx in surfel_idxs {
                Self::absorb(ton, surf.samples[idx].data_mut(), count_weight);
            }

            Self::bounce(surf, tracer, ton, hit, next_motion_type);
        }
    }

    fn bounce(surf: &mut Surface, tracer: &Tracer, ton: &mut Ton, hit: Hit, motion_type: MotionType) {
        let Hit { intersection_point, triangle, incoming_direction } = hit;

        let hit = match motion_type {
            MotionType::Straight => {
                // Regard surface as completely diffuse:
                // Uniformly sample upper Z hemisphere and transform to world space
                // by multiplying with TBN matrix of triangle.
                // This assumes CCW winding order for all triangles since normals, binormals
                // and tangents are calculated from vertices.
                let outgoing_world = triangle.tangent_to_world_matrix() * UnitHemisphere::PosZ.uniform();
                tracer.trace_straight(intersection_point, outgoing_world)
            },
            MotionType::Parabolic => {
                // Also sample diffuse just like in straight
                let outgoing_world = triangle.tangent_to_world_matrix() * UnitHemisphere::PosZ.uniform();
                tracer.trace_parabolic(intersection_point, outgoing_world, ton.parabola_height)
            },
            MotionType::Flow => {
                let up = triangle.normal();
                // TODO project Y unit vector (configurable)
                // REVIEW maybe project_onto_tangential_plane should return an option for the perpendicular case
                let mut flow_direction = triangle.project_onto_tangential_plane(match &ton.flow_direction {
                    &FlowDirection::Incident => incoming_direction,
                    &FlowDirection::Static(global_flow_direction) => global_flow_direction
                });

                // If projected flow direction is zero, it is the result of a failed attempt to project
                // a perpendicular incoming direction, as a fallback, just project a diffusely sampled
                // vector (until a non-zero projection turns up)
                if flow_direction.is_zero() {
                    flow_direction = loop {
                        let fallback_flow_direction = triangle.project_onto_tangential_plane(triangle.uniform());
                        if !fallback_flow_direction.is_zero() {
                            break fallback_flow_direction;
                        }
                    };
                }

                tracer.trace_flow(intersection_point, up, flow_direction, ton.flow_distance)
            },
            MotionType::Settled => unreachable!()
        };

        // REVIEW Some if did not fly out of scene without hitting anything.
        //        In the none case, the material is lost. Is this ok?
        if let Some(hit) = hit {
            Self::interact(surf, tracer, ton, hit);
        }
    }

    fn select_motion_type(ton: &Ton) -> MotionType {
        let mut rng = rand::thread_rng();
        let random : f32 = rng.gen();

        let &Ton { p_straight, p_parabolic, p_flow, .. } = ton;

        if random < p_straight {
            MotionType::Straight
        } else if random < (p_straight + p_parabolic) {
            MotionType::Parabolic
        } else if random < (p_straight + p_parabolic + p_flow) {
            MotionType::Flow
        } else {
            MotionType::Settled
        }
    }

    fn deteriorate(ton: &mut Ton, surfel: &SurfelData) {
        ton.p_straight = (ton.p_straight - surfel.delta_straight).max(0.0);
        ton.p_parabolic = (ton.p_parabolic - surfel.delta_parabolic).max(0.0);
        ton.p_flow = (ton.p_flow + ton.p_parabolic - surfel.delta_flow).max(0.0);

        // REVIEW why this equation for flow? why doesn't it deteriorate like the others
        //        the sum could be larger than 1.0
        //
        // Alternative (sane?) version:
        // ton.p_flow = (ton.p_flow - surfel.delta_flow).max(0.0);

        if (ton.p_straight + ton.p_parabolic + ton.p_flow) > 1.0 {
            warn!("WARN: The flow equation from Chen et. al. yields probability sums > 1.0, fixing it by reducing flow probability");
            warn!("Ton: {:?}", ton);
            warn!("Surfel: {:?}", surfel);
            ton.p_flow -= ton.p_straight + ton.p_parabolic + ton.p_flow - 1.0
        }
    }

    /// Makes the ton pick up material from a surfel it is interacting with.
    /// The pick up rate can also be negative, the ton then deposits material on contact
    /// instead of accumulating.
    fn absorb(ton: &mut Ton, interacting_surfel: &mut SurfelData, count_weight: f32) {
        assert_eq!(
            interacting_surfel.substances.len(), ton.substances.len(),
            "Surfel and ton have unequal amount of materials, cannot transport"
        );

        let material_transports = ton.pickup_rates.iter()
            .zip(
                ton.substances
                    .iter_mut()
                    .zip(
                        interacting_surfel.substances.iter_mut()
                    )
            );

        for (ref pickup_rate, (ref mut ton_material, ref mut surfel_material)) in material_transports {
            // pickup rate gets equally distributed between all interacting surfels
            let pickup_rate = count_weight * *pickup_rate;

            let transport_amount = pickup_rate * if pickup_rate >= 0.0 {
                **surfel_material
            } else {
                **ton_material
            };

            **surfel_material = (**surfel_material - transport_amount).max(0.0);
            **ton_material = (**ton_material + transport_amount).max(0.0);
        }
    }

    /// Deposits the materials in the ton in the interacting surfel, not mutating
    /// the ton
    fn deposit(ton: &Ton, interacting_surfel: &mut SurfelData, count_weight: f32) {
        assert_eq!(
            interacting_surfel.substances.len(), ton.substances.len(),
            "Surfel and ton have unequal amount of materials, cannot transport"
        );

        let material_transports = interacting_surfel.deposition_rates.iter()
            .zip(
                ton.substances
                    .iter()
                    .zip(
                        interacting_surfel.substances.iter_mut()
                    )
            );

        for (ref deposition_rate, (ref ton_material, ref mut surfel_material)) in material_transports {
            // pickup rate gets equally distributed between all interacting surfels
            let deposition_rate = count_weight * *deposition_rate;
            let transport_amount = deposition_rate * **ton_material;
            **surfel_material = (**surfel_material + transport_amount).max(0.0);
        }
    }
}
