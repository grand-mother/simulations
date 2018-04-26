#!/usr/bin/env python
## =============================================================================
#$ -N zhaires-trend50
#$ -P P_trend
#$ -j y
#$ -cwd
#$ -notify
#$ -l ct=48:00:00
#$ -l vmem=3.0G
#$ -l fsize=3.0G
#$ -l sps=1
## =============================================================================
import json
import os
import random
import shutil
import struct
import subprocess
import sys
import time

import numpy as np

# Aires binary.
AIRESBIN = "/afs/in2p3.fr/throng/trend/soft/sim/aires/bin/Aires"

# Output dir for Aires results.
DATADIR = "/sps/hep/trend/zhaires/trend-50/without-freq-sandra-iron"

# Energies to sample, in eV.
#ENERGIES = (8E+16, 8E+16)
#ENERGIES = (5E+16, 1E+17, 3E+17, 5E+17, 1E+18, 3E+18)
ENERGIES = (7E+16,7E+16)

# Range for the zenith angle.
ZENITH_RANGE = (0., 80.)

# Cut on the range.
RANGE_CUT = 600.

# Minimim number of antennas within the range cut.
MIN_ANTENNAS = 5

# Reference altitude for Ulastai.
ZREF = 2650.

# Geomagnetic field (Amplitude [uT], inclination [deg], declination [deg]).
GEOMAGNET = (56.5, 63.18, 2.72)

# Configure the environment.
TAG = os.getenv("JOB_ID")
if TAG is None:
	TAG = ""
	TMPDIR = "/tmp/" + os.getenv("USER")
	if not os.path.exists(TMPDIR):
		os.makedirs(TMPDIR)
else:
	TMPDIR = os.getenv("TMPDIR")
	task = os.getenv("TASK_ID")
	if task:
		TAG += "{:03d}".format(int(task) - 1)
ROOTDIR = os.getcwd()
os.chdir(TMPDIR)

# Task index counter.
TASK_INDEX = 0

# Antennas positions in TREND coordinates system.
ANTENNAS = np.array([
	[101, 1466.7, 3.8, 2653.0],
	[102, 1464.9, 68.2, 2653.0],
	[103, 1594.0, -205.2, 2653.0],
	[104, 1592.5, 194.8, 2656.0],
	[105, 1742.1, -80.0, 2656.0],
	[106, 1742.4, 69.8, 2657.0],
	[107, 1892.7, -202.3, 2659.0],
	[108, 1903.0, 152.5, 2659.0],
	[109, 2226.7, -5.8, 2668.0],
	[110, 2042.6, -81.1, 2663.0],
	[111, 2042.6, 68.3, 2664.0],
	[112, 2193.8, -222.3, 2668.0],
	[113, 2197.6, 188.7, 2667.0],
	[114, 2347.4, -81.1, 2673.0],
	[115, 2348.9, 68.5, 2671.0],
	[116, 2498.2, 192.2, 2674.0],
	[117, 2495.7, -206.3, 2684.0],
	[118, 2507.7, -6.4, 2677.0],
	[119, 2646.1, -82.0, 2687.0],
	[120, 2646.5, 67.8, 2680.0],
	[121, 1319.0, -243.0, 2646.8],
	[122, 1299.7, 187.1, 2646.0],
	[123, 1127.4, 60.0, 2645.0],
	[124, 1133.2, -89.3, 2645.0],
	[125, 983.6, 170.4, 2647.0],
	[126, 989.3, -213.0, 2647.0],
	[127, 848.5, -86.6, 2644.0],
	[128, 848.4, 62.0, 2643.0],
	[129, 699.1, 194.5, 2643.0],
	[130, 699.4, -205.8, 2643.0],
	[131, 549.3, 66.3, 2640.0],
	[132, 548.7, -79.9, 2640.0],
	[133, 479.9, 165.6, 2639.5],
	[134, 462.3, -130.6, 2640.0],
	[135, 369.6, 80.6, 2637.0],
	[136, 304.0, -112.0, 2636.0],
	[137, 232.0, 26.6, 2633.0],
	[138, 224.6, -45.3, 2634.0],
	[140, 100.6, -104.6, 2635.0],
	[148, 49.9, 559.5, 2640.5],
	[149, -19.9, 530.3, 2640.2],
	[150, 90.8, 460.2, 2638.3],
	[151, 23.7, 449.4, 2639.4],
	[152, -33.4, 385.6, 2637.4],
	[153, 82.8, 342.0, 2635.7],
	[154, -33.6, 256.6, 2633.1],
	[155, 58.1, 190.1, 2632.0],
	[156, -20.6, -1, 2632.7],
	[157, 22.5, -100.7, 2634.8],
	[158, 72, -199, 2637.]])

def urandom():
	"""Generate a random float in [0,1] using os.urandom."""
	return struct.unpack("<L", os.urandom(4))[0] / float((1 << 32) - 1)

def generate_input(energy, azimuth, zenith, core, antennas=None):
	"""Generate the input stream for ZHAIRES."""

	# Check the inputs.
	if antennas is None:
		antennas = ANTENNAS

	# Rotation matrix in (x, y) for antennas.
	a = GEOMAGNET[2] * np.pi / 180.
	s, c = np.sin(a), np.cos(a)
	R = np.array(((s, c), (-c, s)))

	# Format the stream.
	stream = [
		"TaskName {:s}{:03d}".format(TAG, TASK_INDEX),
		"PrimaryParticle iron",
		"PrimaryEnergy {:.5E} eV".format(energy),
		"PrimaryZenAngle {:.5f} deg".format(zenith),
		"PrimaryAzimAngle {:.5f} deg Magnetic".format(azimuth),
		"ForceModelName SIBYLL"
	]
	for a in antennas:
		xT = a[1] - core[0]
		yT = a[2] - core[1]
		xZ, yZ = np.dot(R, (xT, yT))
		stream.append("AddAntenna {:.5f} {:.5f} {:.5f}".format(
			xZ, yZ, a[3] - ZREF))
	stream += [
		"TotalShowers 1",
		"RunsPerProcess Infinite",
		"ShowersPerRun 1",
		"RandomSeed {:.12E}".format(urandom()),
		"InjectionAltitude 100 km",
		"Atmosphere 1",
		"AddSite Ulastai 42.55 deg 86.68 deg {:.3f} m".format(ZREF),
		"Site Ulastai",
		"Date 1985 10 26",
		"GeomagneticField On",
		"GeomagneticField {:.4f} uT {:.2f} deg {:.2f} deg".format(
		    *GEOMAGNET),
		"ObservingLevels 510 0 g/cm2   750 g/cm2",
		"PerShowerData Full",
		"SaveNotInFile lgtpcles All",
		"SaveNotInFile grdpcles All",
		"RLimsFile grdpcles 0.000001 m 10 km",
		"ResamplingRatio 100",
		"RLimsTables 10 m 10 km",
		"ELimsTables 2 MeV 1 TeV",
		"ExportTables 5501 Opt a",
		"ExportTables 1293 Opt a",
		"ZHAireS On",
		"FresnelTime On",
		"TimeDomainBin 1 ns",
		"ElectronCutEnergy 3 MeV",
		"ElectronRoughCut 3 MeV",
		"GammaCutEnergy 3 MeV",
		"GammaRoughCut 3 MeV",
		"ThinningEnergy 1.e-4 Relative",
		"ThinningWFactor 0.06"
	]
	return "\n".join(stream)

def sample_shower():
	# Draw the main shower parameters.
	energy = random.choice(ENERGIES)
	deg = np.pi / 180.
	azimuth = random.uniform(0., 360.)
	zenith = np.arccos(random.uniform(np.cos(ZENITH_RANGE[1] * deg),
		np.cos(ZENITH_RANGE[0] * deg))) / deg
	#phi = (90. - azimuth) * deg
	phi=(azimuth+90)*deg #SL
        if phi>=2*np.pi:
            phi=phi-2*np.pi
        theta = zenith * deg
	u = np.array((np.cos(phi) * np.sin(theta),
		np.sin(phi) * np.sin(theta), np.cos(theta)))

	# Sample the core position.
	b = np.mean(ANTENNAS[:,1:], axis=0)
	dz = np.max(np.absolute(ANTENNAS[:,3] - b[2]))
	dx = np.max(np.absolute(ANTENNAS[:,1] - b[0])) + dz + RANGE_CUT
	dy = np.max(np.absolute(ANTENNAS[:,2] - b[1])) + dz + RANGE_CUT
	exposure = dx * dy * 2. * np.pi * (np.cos(ZENITH_RANGE[0]) -
		np.cos(ZENITH_RANGE[1]))
		
	events = 0
	while True:
		events += 1
		x = random.uniform(b[0] - dx, b[0] + dx) / u[2]
		y = random.uniform(b[1] - dy, b[1] + dy) / u[2]
		r0 = np.array((x, y, ZREF))
		# Check the antennas at range.
		antennas = []
		for i, antenna in enumerate(ANTENNAS):
			dr = antenna[1:] - r0
			rho2 = np.linalg.norm(dr)**2 - np.dot(dr, u)**2
			if rho2 <= RANGE_CUT ** 2:
				antennas.append(antenna)
		if len(antennas) >= MIN_ANTENNAS:
			break		
	antennas = np.array(antennas)
	
	# Return the sampled shower and the selected antennas.
	return energy, azimuth, zenith, (x, y), events, exposure, antennas

def log(txt):
	print txt
	sys.stdout.flush()

def run(energy, azimuth, zenith, core, events, exposure, antennas=None):
	"""Run a single ZHAIRES shower + fields computation."""
	global TASK_INDEX

	# Format the inputs for ZHAIRES.
	instream = generate_input(energy, azimuth, zenith, core, antennas)
	
	# Dump a summary.
	outdir = "{:s}{:03d}".format(TAG, TASK_INDEX)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	summary = {"energy": energy, "azimuth": azimuth,
		"zenith": zenith, "core": core, "generated": events,
		"exposure": exposure, "antennas": map(int, antennas[:,0]),
		"input": instream}
	summary = json.dumps(summary, sort_keys=True, indent=4,
		separators=(',', ': '))
	with open(os.path.join(outdir, "summary.json"), "w+") as f:
		f.write(summary)
		
	# Run ZHAIRES.
	p = subprocess.Popen(AIRESBIN, stdin=subprocess.PIPE,
		stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	log("* Processing " + outdir)
	t0 = time.time() 
	stdout, stderr = p.communicate(instream)
	log("  --> Done in {:.0f} s".format(time.time() - t0))

	# Copy back the log and the results.
	etag = "{:.0e}".format(energy).replace("+", "")
	datadir = os.path.join(DATADIR, etag)
	if not os.path.exists(datadir): os.makedirs(datadir)
	with open(os.path.join(datadir, "{:}.log".format(outdir)), "w+") as f:
		f.write(stdout)
	shutil.move("Aires.dirs", outdir)
	shutil.move("Aires.status", outdir)
	shutil.move("timefresnel-root.dat", outdir)
	tgz = outdir + ".tgz"
	cmd = "tar -czf {:} {:}".format(tgz, outdir)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
		stderr=subprocess.PIPE, shell=True)
	_, _ = p.communicate()
	shutil.move(tgz, os.path.join(datadir, tgz))

	# Check for any error.
	if stderr:
		raise RuntimeError(stderr)

	# Update the task index.
	TASK_INDEX += 1

if __name__ == "__main__":
	# Process as many shower + fields as possible.
	while True:
		run(*sample_shower())
