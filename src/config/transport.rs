use self::Transport::*;
use std::default::Default;
use transport::{transport, Absorb, AbsorbThenDeposit, Deposit, Transport as Inner};

/// Specifies when and in which direction substance is transported.
pub enum Transport {
    /// Absorb on bounce, deposit on settle.
    Classic(Inner<Absorb, Deposit>),
    /// On both bounce and settle, absorb first and then deposit.
    Consistent(Inner<AbsorbThenDeposit, AbsorbThenDeposit>),
}

impl Transport {
    pub fn classic() -> Self {
        Classic(transport())
    }

    pub fn consistent() -> Self {
        Consistent(transport())
    }
}

impl Default for Transport {
    fn default() -> Self {
        Self::classic()
    }
}
