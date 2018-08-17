use motion::MotionType;
use std::marker::PhantomData;
use surf::Surfel;
use ton::Ton;
use SurfelData;

/// Performs substance transport for use in surface contacts.
pub struct Transport<B, S> {
    bounce: PhantomData<B>,
    settle: PhantomData<S>,
}

pub fn transport<B: Rule, S: Rule>() -> Transport<B, S> {
    Transport {
        bounce: PhantomData,
        settle: PhantomData,
    }
}

pub fn classic_transport() -> Transport<Absorb, Deposit> {
    transport::<Absorb, Deposit>()
}

pub fn consistent_transport() -> Transport<AbsorbThenDeposit, AbsorbThenDeposit> {
    transport::<AbsorbThenDeposit, AbsorbThenDeposit>()
}

impl<B, S> Transport<B, S>
where
    B: Rule,
    S: Rule,
{
    pub fn perform<V>(
        &self,
        ton: &mut Ton,
        surface: &mut Vec<Surfel<V, SurfelData>>,
        interaction_info: &(MotionType, Vec<usize>),
    ) {
        let (next_motion_type, surfel_idxs) = interaction_info;
        let count_weight = (surfel_idxs.len() as f32).recip();

        match next_motion_type {
            // settle substance exchange
            MotionType::Settled => {
                for idx in surfel_idxs {
                    S::transport(ton, surface[*idx].data_mut(), count_weight);
                }
            }
            // non-settle substance exchange, a.k.a bounce
            _ => {
                for idx in surfel_idxs {
                    B::transport(ton, surface[*idx].data_mut(), count_weight);
                }
            }
        }
    }
}

pub trait Rule {
    fn transport(ton: &mut Ton, interacting_surfel: &mut SurfelData, count_weight: f32);
}

pub struct Absorb;

pub struct Deposit;

impl Rule for Absorb {
    fn transport(ton: &mut Ton, interacting_surfel: &mut SurfelData, count_weight: f32) {
        absorb(ton, interacting_surfel, count_weight);
    }
}

impl Rule for Deposit {
    fn transport(ton: &mut Ton, interacting_surfel: &mut SurfelData, count_weight: f32) {
        deposit(ton, interacting_surfel, count_weight);
    }
}

pub struct AbsorbThenDeposit;
impl Rule for AbsorbThenDeposit {
    fn transport(ton: &mut Ton, interacting_surfel: &mut SurfelData, count_weight: f32) {
        absorb(ton, interacting_surfel, count_weight);
        deposit(ton, interacting_surfel, count_weight);
    }
}

/// Deposits the materials in the ton in the interacting surfel, not mutating
/// the ton
fn deposit(ton: &Ton, interacting_surfel: &mut SurfelData, count_weight: f32) {
    assert_eq!(
        interacting_surfel.substances.len(),
        ton.substances.len(),
        "Surfel and ton have unequal amount of materials, cannot transport"
    );

    let material_transports = interacting_surfel.deposition_rates.iter().zip(
        ton.substances
            .iter()
            .zip(interacting_surfel.substances.iter_mut()),
    );

    for (ref deposition_rate, (ref ton_material, ref mut surfel_material)) in material_transports {
        // pickup rate gets equally distributed between all interacting surfels
        let deposition_rate = count_weight * *deposition_rate;
        let transport_amount = deposition_rate * **ton_material;
        **surfel_material = (**surfel_material + transport_amount).max(0.0);
    }
}

/// Makes the ton pick up material from a surfel it is interacting with.
/// The pick up rate can also be negative, the ton then deposits material on contact
/// instead of accumulating.
fn absorb(ton: &mut Ton, interacting_surfel: &mut SurfelData, count_weight: f32) {
    assert_eq!(
        interacting_surfel.substances.len(),
        ton.substances.len(),
        "Surfel and ton have unequal amount of materials, cannot transport"
    );

    let material_transports = ton.pickup_rates.iter().zip(
        ton.substances
            .iter_mut()
            .zip(interacting_surfel.substances.iter_mut()),
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
