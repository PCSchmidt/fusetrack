//! Output stage: publish a confirmed track to a Cursor on Target (CoT) event
//! for TAK, via the [`uas-cot`](https://docs.rs/uas-cot) crate.
//!
//! The tracker produces confirmed tracks in the local Cartesian frame
//! (metres east and north of the sensor). A common operating picture needs
//! them as WGS84 CoT events. [`track_to_cot`] converts one confirmed track
//! into an unknown-UAS CoT event (`a-u-A-M-F-Q`) that ATAK, WinTAK and a TAK
//! Server accept.

use uas_cot::{Event, Point, Track as CotTrack};

use crate::tracker::{Track, TrackStatus};

/// WGS84 origin of the tracker's local Cartesian frame: the point that
/// `[x, y] = [0, 0]` maps to, with `+x` east and `+y` north.
#[derive(Clone, Copy, Debug)]
pub struct Origin {
    /// Latitude of the frame origin, degrees.
    pub lat: f64,
    /// Longitude of the frame origin, degrees.
    pub lon: f64,
    /// Height of the local plane above the WGS84 ellipsoid, metres.
    pub hae: f64,
}

const WGS84_A: f64 = 6_378_137.0;
const WGS84_E2: f64 = 6.694_379_990_14e-3;

impl Origin {
    /// WGS84 latitude/longitude (degrees) of a local `[east, north]` offset
    /// in metres, on the tangent plane at this origin.
    fn to_wgs84(&self, east: f64, north: f64) -> (f64, f64) {
        let lat = self.lat.to_radians();
        let s = lat.sin();
        let w = (1.0 - WGS84_E2 * s * s).sqrt();
        // Meridional and prime-vertical radii of curvature at the origin.
        let m = WGS84_A * (1.0 - WGS84_E2) / (w * w * w);
        let n = WGS84_A / w;
        let dlat = (north / m).to_degrees();
        let dlon = (east / (n * lat.cos())).to_degrees();
        (self.lat + dlat, self.lon + dlon)
    }
}

/// Turn a confirmed track into an unknown-UAS CoT event (`a-u-A-M-F-Q`) at
/// `time_us` (microseconds since 1970 UTC), stale after `stale_secs`.
///
/// Returns `None` for a tentative track: only confirmed tracks are published
/// to the common operating picture. Course (degrees true) and speed (m/s)
/// come from the track's velocity.
pub fn track_to_cot(
    track: &Track,
    origin: &Origin,
    time_us: u64,
    stale_secs: u64,
) -> Option<Event> {
    if track.status != TrackStatus::Confirmed {
        return None;
    }
    let x = track.filter.state();
    let (lat, lon) = origin.to_wgs84(x[0], x[1]);
    let mut e = Event::new(
        format!("fusetrack-{}", track.id),
        "a-u-A-M-F-Q",
        time_us,
        stale_secs,
    );
    e.point = Point::new(lat, lon, origin.hae);
    let (vx, vy) = (x[2], x[3]);
    let speed = vx.hypot(vy);
    if speed > 0.0 {
        let course = vx.atan2(vy).to_degrees().rem_euclid(360.0);
        e.track = Some(CotTrack { course, speed });
    }
    Some(e)
}
