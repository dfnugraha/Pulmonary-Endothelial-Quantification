// PLVAP: 75% DAPI/ERG gate + 5% local PLVAP/P99 cell counts; THREE final COLORED QC images
// Adapted for OlyVIA 3-plane 8/16-bit raw TIFF, including 1539 x 1539 images
// Derived from the preceding full-resolution PLVAP cell-counting macro; requires Fiji validation
// Reads physical pixel size from TIFF calibration. Does not assume a 512 px field.
//
// Logic:
// DAPI nuclei -> ERG-positive nuclei -> expand by distance -> PLVAP mask -> accepted PLVAP = PLVAP_mask AND ERG_near_zone
//
// Expected raw Olympus VSI channel order for Sophie's PLVAP test images:
// C1 = DAPI
// C2 = Alexa488 / ERG
// C3 = Alexa594 (AF568 acquisition project) / PLVAP
//
// IMPORTANT: display pseudocolors in OlyVIA do not change this raw channel order.
//
// Main fixes:
// Final association now checks accepted PLVAP pixels within a fixed
// neighborhood of each original DAPI+/ERG+ nuclear ROI, without using
// ROI Enlarge (which failed for ROI index 2).
// The previous version left the last nucleus ROI selected.
// In Fiji, Duplicate can duplicate only the active selection, not the full image.
// This made ERG_near_zone wrong. This version clears selections before duplicating masks.

requires("1.53");

orig = getTitle();
prefix = "ERGdistPLVAPstrict_assocERG_" + d2s(getTime(), 0);

// ================= USER SETTINGS =================

cERG   = 2;
cDAPI  = 1;
cPLVAP = 3;

// Pixel size will be read from ImageJ calibration in the TIFF.

// Balanced distance around ERG-positive nuclei
nearDistance_um = 5;

// Background subtraction
// Spatially scaled from full-res 30/50/100 px at 0.325012 um/px
rollingBallDAPI  = 10;
rollingBallERG   = 17;
rollingBallPLVAP = 33;

// Spatially scaled Gaussian sigmas (full-res: 1/1/0.3 px)
gaussianSigmaDAPI  = 0.33;
gaussianSigmaERG   = 0.33;
gaussianSigmaPLVAP = 0.10;

// ERG threshold carried forward from the supplied 512 x 512 macro.
// The 488/ERG channel contains substantial diffuse septal signal.
// Strict test: mean + 2.5 SD after background subtraction, and a DAPI nucleus
// must have >=75% overlap with above-threshold ERG pixels. This is not a
// biological truth: inspect nuclear QC masks, and compare with previous settings.
// Both criteria are configurable below.
dapiThresholdMethod  = "Otsu";
ergThresholdMethod   = "Manual_SD";
ergSDmultiplier         = 2.5;      // stricter ERG threshold, previously 2 SD
plvapThresholdMethod = "Manual_SD";  // stricter: mean + SD multiplier threshold

// DAPI nucleus filtering, um^2
dapiMinArea = 15;
dapiMaxArea = 300;

// Strict ERG-positive nucleus definition
ergNuclearFractionThreshold = 0.75;  // >=75% of DAPI nucleus must exceed ERG intensity threshold

// PLVAP object filtering, um^2
plvapMinArea = 0.20;   // small puncta allowed
plvapMaxArea = 25;     // strict: excludes broad parenchymal/septal regions

// PLVAP puncta/stringency settings
// Increase to 3.5 or 4.0 if parenchyma is still counted.
// Decrease to 2.0 if true PLVAP puncta disappear.
plvapSDmultiplier = 3.0;

// PLVAP-positive cell definition matches the CAR4 cell-counting macro.
// This is independent of the PLVAP particle (mean+3SD) threshold.
plvapCellMinLocalFraction = 0.05;
plvapCellBgExclusion_um = 10;

// Circularity filter for PLVAP objects.
// 0.00 allows long structures. 0.10-0.20 removes many long septal lines.
// If true PLVAP is elongated, use 0.00.
plvapCircularityMin = 0.10;

dapiWatershed = true;
// Keep Despeckle off to avoid removing small PLVAP structures.
plvapDoDespeckle = false;

// ==================================================

// ==================================================
// TIFF input normalization and spatial calibration
// Exported OlyVIA TIFF may open as a 3-page Z-stack, not 3 C channels.
// Convert the pages to a 3-channel hyperstack while preserving raw values.
// Also accept an image that already opens as C=3, Z=1, T=1.
// ==================================================
getDimensions(imgWidth, imgHeight, imgChannels, imgSlices, imgFrames);
imageBits = bitDepth();

if (imageBits != 8 && imageBits != 16)
    exit("Input must be 8-bit or 16-bit grayscale TIFF planes; found " + imageBits + "-bit. Do not use an RGB-merged image.");

if (imgChannels == 1 && imgSlices == 3 && imgFrames == 1) {
    run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames=1 display=Grayscale");
    orig = getTitle();
    getDimensions(imgWidth, imgHeight, imgChannels, imgSlices, imgFrames);
    print("Converted three TIFF pages to C1/C2/C3 hyperstack.");
}

if (imgChannels != 3 || imgSlices != 1 || imgFrames != 1)
    exit("Open the original 3-plane raw TIFF stack (C=3 or Z=3). Current dimensions: C=" + imgChannels +
         ", Z=" + imgSlices + ", T=" + imgFrames + ". Do not use the split channel windows.");

getPixelSize(scaleUnit, pixelWidth_um, pixelHeight_um);
scaleUnitLower = toLowerCase(scaleUnit);
if (scaleUnitLower != "um" && scaleUnitLower != "µm" && scaleUnitLower != "micron" &&
    scaleUnitLower != "microns" && scaleUnitLower != "micrometer" &&
    scaleUnitLower != "micrometers" && scaleUnitLower != "micrometre" &&
    scaleUnitLower != "micrometres")
    exit("Spatial calibration is missing or not in micrometers (unit=" + scaleUnit + "). Set it from trusted acquisition metadata using Image > Properties before running.");

if (pixelWidth_um <= 0 || pixelHeight_um <= 0 ||
    abs(pixelWidth_um - pixelHeight_um) > 0.01 * pixelWidth_um)
    exit("Unexpected/non-square micrometer calibration: X=" + pixelWidth_um + ", Y=" + pixelHeight_um + ". Verify the TIFF metadata.");

pixelSize_um = pixelWidth_um;
fieldWidth_um = imgWidth * pixelWidth_um;
fieldHeight_um = imgHeight * pixelHeight_um;

// Retain the ORIGINAL physical preprocessing radii of the 512-pixel macros.
// The new 1539-pixel export reports ~0.977 µm/pixel; it covers ~1503.6 µm,
// so tripling each radius merely because pixel dimensions tripled would be wrong.
referencePixelSize_um = 500.2 / 512;
rollingBallDAPI  = round(rollingBallDAPI * referencePixelSize_um / pixelSize_um);
rollingBallERG   = round(rollingBallERG * referencePixelSize_um / pixelSize_um);
// Marker-specific rolling balls / Gaussian blur parameters are set below.
if (rollingBallDAPI < 1) rollingBallDAPI = 1;
if (rollingBallERG < 1) rollingBallERG = 1;
rollingBallPLVAP = round(rollingBallPLVAP * referencePixelSize_um / pixelSize_um);
if (rollingBallPLVAP < 1) rollingBallPLVAP = 1;
gaussianSigmaDAPI = gaussianSigmaDAPI * referencePixelSize_um / pixelSize_um;
gaussianSigmaERG = gaussianSigmaERG * referencePixelSize_um / pixelSize_um;
gaussianSigmaPLVAP = gaussianSigmaPLVAP * referencePixelSize_um / pixelSize_um;

run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
run("Set Measurements...", "area mean min max integrated redirect=None decimal=3");
setOption("BlackBackground", true);
setBackgroundColor(0, 0, 0);
setForegroundColor(255, 255, 255);

ergDilatePixels = round(nearDistance_um / pixelSize_um);

print("\\Clear");
print("Image: " + orig);
print("Unique run prefix: " + prefix);
print("Image size: " + imgWidth + " x " + imgHeight);
print("Channels: " + imgChannels + ", Slices: " + imgSlices + ", Frames: " + imgFrames);
print("Input: " + imgWidth + " x " + imgHeight + ", C=" + imgChannels + ", " + imageBits + "-bit");
print("Calibrated field: " + fieldWidth_um + " x " + fieldHeight_um + " um");
print("Pixel size: " + pixelSize_um + " um/pixel");
print("Raw channel order: C1 DAPI / C2 ERG / C3 PLVAP (VERIFY!)");
print("ERG expansion distance: " + nearDistance_um + " um");
print("ERG expansion in pixels: " + ergDilatePixels);
print("ERG association uses direct pixel-neighborhood tests; no ROI enlargement.");

// Window names
ergRaw = prefix + "_ERG_raw";
dapiRaw = prefix + "_DAPI_raw";
plvapRaw = prefix + "_PLVAP_raw";
dapiMask = prefix + "_DAPI_mask";
ergThresh = prefix + "_ERG_thresh_mask";
ergDapiOverlapMask = prefix + "_ERG_above_threshold_INSIDE_DAPI_QC";
ergTwoColorQC = prefix + "_DAPI_blue_ERG_green_raw_QC";
ergNuclearTwoColorQC = prefix + "_DAPI_blue_ERG_nuclear_red_QC";
ergPosMask = prefix + "_ERGpos_nuclei_mask";
ergZone = prefix + "_ERG_near_zone";
plvapBgsub = prefix + "_PLVAP_bgsub";
plvapMask = prefix + "_PLVAP_threshold_mask";
plvapFilteredMask = prefix + "_PLVAP_candidate_filtered_mask";
acceptedPreMask = prefix + "_Accepted_PLVAP_pre_filter";
acceptedMask = prefix + "_Accepted_PLVAP_mask";
acceptedRawMap = prefix + "_Accepted_PLVAP_raw_intensity_map";
acceptedBgsubMap = prefix + "_Accepted_PLVAP_bgsub_intensity_map";
plvapAssocERGMask = prefix + "_PLVAP_associated_ERGpos_nuclei_mask";
plvapNearForERG = prefix + "_PLVAP_near_for_ERG_count";
associationOverlapQC = prefix + "_PLVAPnear_AND_ERGpos_QC";
tripleCandidateQC = prefix + "_DAPI_blue_ERG_nuclear_red_PLVAP_candidates_green_QC";
tripleAcceptedQC = prefix + "_DAPI_blue_ERGpos_red_ACCEPTED_PLVAP_green_QC";

// ==================================================
// 1. Duplicate raw channels
// ==================================================

selectWindow(orig);
run("Select None");
Stack.setChannel(cERG);
run("Duplicate...", "title=" + ergRaw + " duplicate channels=" + cERG);

selectWindow(orig);
run("Select None");
Stack.setChannel(cDAPI);
run("Duplicate...", "title=" + dapiRaw + " duplicate channels=" + cDAPI);

selectWindow(orig);
run("Select None");
Stack.setChannel(cPLVAP);
run("Duplicate...", "title=" + plvapRaw + " duplicate channels=" + cPLVAP);

// ==================================================
// 2. DAPI nuclei mask and ROIs
// ==================================================

selectWindow(dapiRaw);
run("Select None");
run("Duplicate...", "title=" + dapiMask);

selectWindow(dapiMask);
run("Subtract Background...", "rolling=" + rollingBallDAPI);
run("Gaussian Blur...", "sigma=" + gaussianSigmaDAPI);
setAutoThreshold(dapiThresholdMethod + " dark");
run("Convert to Mask");

getStatistics(areaDapiCheck, meanDapiCheck);
if (meanDapiCheck > 127) {
    run("Invert");
}

if (dapiWatershed == true) {
    run("Watershed");
}

roiManager("Reset");
run("Analyze Particles...",
    "size=" + dapiMinArea + "-" + dapiMaxArea +
    " circularity=0.10-1.00 show=Nothing exclude add");

nDAPI = roiManager("count");

if (nDAPI == 0) {
    exit("No DAPI nuclei detected. Adjust DAPI threshold or size settings.");
}

print("Total DAPI nuclei detected: " + nDAPI);

// Preserve ALL DAPI nuclear ROIs, including ERG-negative nuclei, for DAPI+/PLVAP+ counts.
dapiRoiZip = getDirectory("temp") + prefix + "_ALL_DAPI_nuclear_ROIs.zip";
roiManager("Save", dapiRoiZip);

// QC1 counts only the nuclear ROIs that passed DAPI area/circularity filtering.
// The original dapiMask can contain thresholded pixels that were not counted.
dapiDetectedQCMask = prefix + "_DAPI_COUNTED_NUCLEI_mask";
newImage(dapiDetectedQCMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
for (di=0; di<nDAPI; di++) {
    selectWindow(dapiMask);
    roiManager("Select", di);
    if (selectionType() == -1) exit("DAPI QC ERROR: missing counted nucleus ROI " + di);
    Roi.getContainedPoints(dapiQX, dapiQY);
    selectWindow(dapiDetectedQCMask);
    for (dj=0; dj<dapiQX.length; dj++) setPixel(dapiQX[dj],dapiQY[dj],255);
}
selectWindow(dapiDetectedQCMask);
run("Select None");

// ==================================================
// 3. ERG threshold mask
// ==================================================

selectWindow(ergRaw);
run("Select None");
run("Duplicate...", "title=" + ergThresh);

selectWindow(ergThresh);
run("Subtract Background...", "rolling=" + rollingBallERG);
run("Gaussian Blur...", "sigma=" + gaussianSigmaERG);

// ERG threshold retained from the original macro. Verify against QC masks.
// The previous Yen threshold included broad septal/background signal.
// Here we retain only pixels above mean + 2 SD of the background-subtracted ERG image.
getStatistics(areaERGbgsub, meanERGbgsub, minERGbgsub, maxERGbgsub, sdERGbgsub);
ergLowerThreshold = meanERGbgsub + ergSDmultiplier * sdERGbgsub;
setThreshold(ergLowerThreshold, maxERGbgsub);
print("ERG manual threshold lower value: " + ergLowerThreshold);
print("ERG threshold rule: mean + " + ergSDmultiplier + " SD");
run("Convert to Mask");

getStatistics(areaErgCheck, meanErgCheck);
if (meanErgCheck > 127) {
    run("Invert");
}

getStatistics(areaErgMask, meanErgMask);
ergMaskFraction = meanErgMask / 255;
print("ERG threshold-mask positive fraction of image: " + ergMaskFraction);

// Double-positive QC: show ONLY above-threshold ERG pixels inside DAPI mask.
// Whole ERG threshold masks can include septal structures outside nuclei.
selectWindow(dapiMask);
run("Select None");
selectWindow(ergThresh);
run("Select None");
imageCalculator("AND create", dapiMask, ergThresh);
rename(ergDapiOverlapMask);
run("Select None");

// ==================================================
// 4. ERG-positive nuclei based on ERG signal inside DAPI ROIs
// ==================================================

newImage(ergPosMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");

nERGpos = 0;
ergPositiveFlag = newArray(nDAPI);


for (i = 0; i < nDAPI; i++) {
    selectWindow(ergThresh);
    roiManager("Select", i);
    getStatistics(areaERGnuc, meanERGnuc);
    ergFractionInNucleus = meanERGnuc / 255;

    if (ergFractionInNucleus >= ergNuclearFractionThreshold) {
        selectWindow(ergPosMask);
        roiManager("Select", i);
        setForegroundColor(255, 255, 255);
        run("Fill");
        ergPositiveFlag[i] = 1;
        nERGpos++;
    }
}

ergPosFraction = nERGpos / nDAPI;

selectWindow(ergPosMask);
run("Select None");

print("ERG-positive nuclei detected: " + nERGpos);
print("ERG-positive fraction of DAPI nuclei: " + ergPosFraction);

if (nERGpos == 0) {
    print("WARNING: No ERG-positive nuclei detected. Lower ergNuclearFractionThreshold or use Otsu.");
}

// NOTE: the previous >80% ERG-positive-nuclei hard-stop has been removed.
// The macro now continues quantification and reports a warning only.
if (ergPosFraction > 0.80) {
    print("WARNING: More than 80% of DAPI nuclei are classified as ERG-positive.");
    print("Quantification will continue. Inspect ERG_thresh_mask and ERGpos_nuclei_mask for QC.");
}

// Preserve exact ERG-positive DAPI ROIs before the ROI Manager is reused
// for PLVAP particle filtering. This lets us later count ERG+ nuclei that
// have accepted PLVAP signal within the same 5-um association distance.
ergRoiZip = getDirectory("temp") + prefix + "_ERGpos_ROIs.zip";
if (nERGpos > 0) {
    for (i=nDAPI-1; i>=0; i--) {
        if (ergPositiveFlag[i] == 0) {
            roiManager("Select", i);
            roiManager("Delete");
        }
    }
    if (roiManager("count") != nERGpos)
        exit("Internal ROI QC failed while preserving ERG-positive nuclei.");
    roiManager("Save", ergRoiZip);
}

// ==================================================
// 4B. Make inspectable two-color QC composites (without changing raw data)
//  DAPI = blue, ERG = green; separate QC shows nuclear ERG overlap in red.
//  The counting uses the original grayscale intensity, NOT the display LUT.
// ==================================================
selectWindow(dapiRaw);
run("Select None");
selectWindow(ergRaw);
run("Select None");
run("Merge Channels...", "c2=[" + ergRaw + "] c3=[" + dapiRaw + "] create keep");
rename(ergTwoColorQC);
run("Select None");

selectWindow(dapiMask);
run("Select None");
selectWindow(ergDapiOverlapMask);
run("Select None");
run("Merge Channels...", "c1=[" + ergDapiOverlapMask + "] c3=[" + dapiMask + "] create keep");
rename(ergNuclearTwoColorQC);
run("Select None");

print("ERG two-color QC (blue=DAPI, green=raw ERG): " + ergTwoColorQC);
print("ERG nuclear overlap QC (blue=DAPI, red=thresholded ERG within DAPI): " + ergNuclearTwoColorQC);
print("Strict double-positive gate: ERG threshold mean + " + ergSDmultiplier + " SD AND >= " + (100 * ergNuclearFractionThreshold) + "% ERG-positive DAPI nuclear area.");

// ==================================================
// 5. Expand ERG-positive nuclei
// ==================================================

selectWindow(ergPosMask);
run("Select None");   // CRITICAL FIX
run("Duplicate...", "title=" + ergZone);

selectWindow(ergZone);
run("Select None");   // CRITICAL FIX before dilating the whole mask

for (i = 0; i < ergDilatePixels; i++) {
    run("Dilate");
}

run("Select None");
getStatistics(areaZoneCheck, meanZoneCheck);
ergZoneFraction = meanZoneCheck / 255;
print("ERG_near_zone positive fraction of image: " + ergZoneFraction);

if (ergZoneFraction > 0.80) {
    print("WARNING: ERG_near_zone covers more than 80% of the image.");
    print("Reduce ERG-positive nuclei count or reduce nearDistance_um.");
}

// ==================================================
// 6. PLVAP mask
// ==================================================

selectWindow(plvapRaw);
run("Select None");
run("Duplicate...", "title=" + plvapBgsub);

selectWindow(plvapBgsub);
run("Subtract Background...", "rolling=" + rollingBallPLVAP);

selectWindow(plvapBgsub);
run("Select None");
run("Duplicate...", "title=" + plvapMask);

selectWindow(plvapMask);

// For puncta-like PLVAP, avoid smoothing too much.
// Strong smoothing can connect diffuse parenchymal background.
run("Gaussian Blur...", "sigma=" + gaussianSigmaPLVAP);

// Manual threshold based on background-subtracted PLVAP intensity.
// This helps avoid counting the whole parenchyma as PLVAP.
getStatistics(areaPLVAPstat, meanPLVAPstat, minPLVAPstat, maxPLVAPstat, stdPLVAPstat);
plvapLowerThreshold = meanPLVAPstat + plvapSDmultiplier * stdPLVAPstat;
setThreshold(plvapLowerThreshold, maxPLVAPstat);
print("PLVAP manual threshold lower value: " + plvapLowerThreshold);
print("PLVAP bgsub mean: " + meanPLVAPstat + ", SD: " + stdPLVAPstat + ", max: " + maxPLVAPstat);
run("Convert to Mask");

// Make sure PLVAP signal is white and background is black.
getStatistics(areaPLVAPcheck, meanPLVAPcheck);
if (meanPLVAPcheck > 127) {
    run("Invert");
}

if (plvapDoDespeckle == true) {
    run("Despeckle");
}

run("Select None");
roiManager("Reset");
// Ask ImageJ to make the mask directly from particles passing the filters.
// Single-pixel particles can have ROI types for which run("Fill") fails.
// The Analyze Particles mask includes these particles without ROI filling.
run("Analyze Particles...",
    "size=" + plvapMinArea + "-" + plvapMaxArea +
    " circularity=" + plvapCircularityMin + "-1.00 show=Masks exclude add");

nPLVAPcandidate = roiManager("count");
if (getTitle() == plvapMask)
    exit("Could not create PLVAP particle mask from Analyze Particles. Verify Fiji supports show=Masks.");
rename(plvapFilteredMask);
run("Select None");

selectWindow(plvapFilteredMask);
run("Select None");
getStatistics(areaPLVAPmask, meanPLVAPmask);
plvapMaskFraction = meanPLVAPmask / 255;

print("Total PLVAP candidate objects before ERG-distance filtering: " + nPLVAPcandidate);
print("Filtered PLVAP candidate-mask positive fraction of image: " + plvapMaskFraction);

// ==================================================
// 7. Accepted PLVAP = FILTERED PLVAP candidates AND ERG_near_zone
// ==================================================

selectWindow(plvapFilteredMask);
run("Select None");
selectWindow(ergZone);
run("Select None");

imageCalculator("AND create", plvapFilteredMask, ergZone);
rename(acceptedPreMask);

selectWindow(acceptedPreMask);
run("Select None");

getStatistics(areaAcceptedPre, meanAcceptedPre);
acceptedPreFraction = meanAcceptedPre / 255;
print("Accepted PLVAP mask positive fraction before final particle filtering: " + acceptedPreFraction);

if (plvapDoDespeckle == true) {
    run("Despeckle");
}

run("Select None");
roiManager("Reset");
// Use the particle-mask output rather than ROI-by-ROI filling (which
// can fail for single-pixel accepted objects).
run("Analyze Particles...",
    "size=" + plvapMinArea + "-" + plvapMaxArea +
    " circularity=" + plvapCircularityMin + "-1.00 show=Masks exclude add");

nAccepted = roiManager("count");
if (getTitle() == acceptedPreMask)
    exit("Could not create accepted PLVAP particle mask from Analyze Particles.");
rename(acceptedMask);
run("Select None");

selectWindow(acceptedMask);
run("Select None");
getStatistics(areaAcceptedMask, meanAcceptedMask);
acceptedMaskFraction = meanAcceptedMask / 255;

print("ERG-distance-associated PLVAP objects: " + nAccepted);
print("Final accepted PLVAP mask positive fraction: " + acceptedMaskFraction);

// ==================================================
// 7B. Create intensity-map QC images
// These show raw/bg-subtracted PLVAP intensity only inside accepted PLVAP objects.
// Everything outside Accepted_PLVAP_mask is cleared to black.
// ==================================================

// If there are no accepted objects, a zero-mask has no area selection.
// Produce empty intensity maps rather than running Create Selection on it.
selectWindow(plvapRaw);
run("Select None");
run("Duplicate...", "title=" + acceptedRawMap);
selectWindow(plvapBgsub);
run("Select None");
run("Duplicate...", "title=" + acceptedBgsubMap);

if (nAccepted > 0) {
    selectWindow(acceptedMask);
    run("Select None");
    run("Create Selection");
    selectWindow(acceptedRawMap);
    run("Restore Selection");
    run("Clear Outside");
    run("Select None");

    selectWindow(acceptedMask);
    run("Select None");
    run("Create Selection");
    selectWindow(acceptedBgsubMap);
    run("Restore Selection");
    run("Clear Outside");
    run("Select None");
} else {
    selectWindow(acceptedRawMap);
    run("Select All");
    run("Clear");
    run("Select None");
    selectWindow(acceptedBgsubMap);
    run("Select All");
    run("Clear");
    run("Select None");
    print("No accepted PLVAP objects: intensity QC maps set to zero.");
}

print("Created intensity QC maps:");
print("1. " + acceptedRawMap);
print("2. " + acceptedBgsubMap);

if (ergZoneFraction > 0.95 && plvapMaskFraction > 0.001 && nAccepted == 0) {
    print("TECHNICAL WARNING: ERG_near_zone is almost full image, PLVAP mask is positive, but accepted PLVAP is 0.");
    print("This should not happen. Close old mask windows and rerun, or check mask polarity.");
}

// ==================================================
// 8. Measure accepted PLVAP intensity
// ==================================================

run("Clear Results");

totalArea = 0;
totalRawIntegratedDensity = 0;
totalBgSubIntegratedDensity = 0;

for (i = 0; i < nAccepted; i++) {
    selectWindow(plvapRaw);
    roiManager("Select", i);
    getStatistics(areaRaw, meanRaw, minRaw, maxRaw, stdRaw);
    rawIntegratedDensity = areaRaw * meanRaw;

    selectWindow(plvapBgsub);
    roiManager("Select", i);
    getStatistics(areaBg, meanBg, minBg, maxBg, stdBg);
    bgSubIntegratedDensity = areaBg * meanBg;

    setResult("Accepted_PLVAP_object_ID", i, i + 1);
    setResult("PLVAP_area_um2", i, areaRaw);
    setResult("Mean_PLVAP_raw", i, meanRaw);
    setResult("Integrated_density_raw", i, rawIntegratedDensity);
    setResult("Mean_PLVAP_bgsub", i, meanBg);
    setResult("Integrated_density_bgsub", i, bgSubIntegratedDensity);

    totalArea += areaRaw;
    totalRawIntegratedDensity += rawIntegratedDensity;
    totalBgSubIntegratedDensity += bgSubIntegratedDensity;
}

updateResults();

// ==================================================
// 8B. PLVAP-associated DAPI+/ERG+ nuclei -- NO ROI enlargement
// ==================================================
// A nucleus is associated when at least one final accepted PLVAP-mask
// pixel is within ergDilatePixels of any pixel inside that original
// DAPI+/ERG+ nuclear ROI. The square pixel neighborhood matches the
// ImageJ binary Dilate operation used to construct ERG_near_zone.
//
// No Enlarge, Restore Selection, or per-particle ROI dilation is used.
// The stored ROIs are used only to find nuclear pixels and to draw the
// positive nucleus into the final QC mask. If there is a mismatch,
// report a technical failure rather than treating it as zero biology.

nPLVAPassocERG = 0;
plvapAssocERGFraction = 0;
associationOverlapQC = prefix + "_Accepted_PLVAP_AND_ERG_near_zone_QC";

newImage(plvapAssocERGMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");

// QC: the accepted particle mask must actually overlap the nuclear near-zone.
selectWindow(acceptedMask);
run("Select None");
getStatistics(aAcceptedQC, mAcceptedQC, minAcceptedQC, maxAcceptedQC);
print("Final PLVAP mask pixel range: " + minAcceptedQC + " to " + maxAcceptedQC);
if (nAccepted > 0 && (maxAcceptedQC == 0 || minAcceptedQC != 0))
    exit("QC FAILURE: positive PLVAP objects were reported, but the final mask is missing/invalid.");

selectWindow(ergZone);
run("Select None");
selectWindow(acceptedMask);
run("Select None");
imageCalculator("AND create", acceptedMask, ergZone);
rename(associationOverlapQC);
run("Select None");
getStatistics(aAssocQC, mAssocQC, minAssocQC, maxAssocQC);
print("Accepted PLVAP overlapping ERG near-zone (mask fraction): " + (mAssocQC / 255));
if (nAccepted > 0 && maxAssocQC == 0)
    exit("QC FAILURE: accepted PLVAP particles do not overlap the ERG near-zone. Inspect mask polarity and image alignment.");

if (nERGpos > 0) {
    roiManager("Reset");
    roiManager("Open", ergRoiZip);
    nERGloaded = roiManager("count");
    if (nERGloaded != nERGpos)
        exit("ROI QC FAILURE: saved ERG-positive nucleus ROIs changed in number.");

    selectWindow(acceptedMask);
    for (i=0; i<nERGloaded; i++) {
        // Read original nucleus pixels, WITHOUT enlarging the ROI.
        selectWindow(acceptedMask);
        roiManager("Select", i);
        if (selectionType() == -1)
            exit("ROI QC FAILURE: original ERG nucleus selection is missing (index " + i + ").");
        Roi.getContainedPoints(nucleusX, nucleusY);
        if (nucleusX.length == 0)
            exit("ROI QC FAILURE: original ERG nucleus has no pixels (index " + i + ").");

        associated = 0;
        for (j=0; j<nucleusX.length && associated==0; j++) {
            nx = nucleusX[j];
            ny = nucleusY[j];
            for (dy=-ergDilatePixels; dy<=ergDilatePixels && associated==0; dy++) {
                py = ny + dy;
                if (py>=0 && py<imgHeight) {
                    for (dx=-ergDilatePixels; dx<=ergDilatePixels && associated==0; dx++) {
                        px = nx + dx;
                        if (px>=0 && px<imgWidth && getPixel(px, py)>0)
                            associated = 1;
                    }
                }
            }
        }

        if (associated==1) {
            nPLVAPassocERG++;
            // Draw from the pixel coordinates retrieved from the ORIGINAL nucleus ROI.
            // Do not transfer an ROI selection between windows: some Fiji ROI types
            // lose the selection when switching images, causing run("Fill") to fail.
            selectWindow(plvapAssocERGMask);
            for (k=0; k<nucleusX.length; k++)
                setPixel(nucleusX[k], nucleusY[k], 255);
        }
    }
    selectWindow(plvapAssocERGMask);
    run("Select None");
    plvapAssocERGFraction = nPLVAPassocERG / nERGpos;
    print("PLVAP-associated DAPI+/ERG+ nuclei: " + nPLVAPassocERG + " / " + nERGpos + " = " + plvapAssocERGFraction);
    print("Association searched original nuclear pixels with a " + ergDilatePixels + "-pixel radius, matching ERG-near-zone dilation.");

    if (nAccepted > 0 && nPLVAPassocERG==0)
        exit("ASSOCIATION QC FAILURE: PLVAP objects overlap the ERG-near-zone, but no original nucleus has a nearby PLVAP pixel. Inspect the nucleus mask/ROIs and binary masks; do not report zero.");
    if (nAccepted==0 && nPLVAPassocERG>0)
        exit("ASSOCIATION QC FAILURE: detected associated nuclei with no accepted PLVAP objects.");

    // Saved ERG-positive ROIs are also needed for the cell-level PLVAP gate.
}

// ==================================================
// 8C. Three-color QC overlays (VISUALIZATION ONLY)
// ==================================================
// These are display composites of the ALREADY CALCULATED masks.
// No thresholds, counts, ROIs, association distance or measurements change.
// Fiji's Merge Channels mapping is: c1=RED, c2=GREEN, c3=BLUE.
// DAPI blue + ERG red overlap looks magenta on screen.
//
// Candidate QC: Shows threshold-positive ERG pixels within DAPI nuclei
// in red and all particle-filtered PLVAP candidates in green.
// Note that red pixels can occur in nuclei NOT passing the 75% ERG gate.
selectWindow(ergDapiOverlapMask);
run("Select None");
selectWindow(plvapFilteredMask);
run("Select None");
selectWindow(dapiMask);
run("Select None");
run("Merge Channels...", "c1=[" + ergDapiOverlapMask + "] c2=[" + plvapFilteredMask + "] c3=[" + dapiMask + "] create keep");
rename(tripleCandidateQC);
run("Select None");

// Accepted QC: Marks complete DAPI+/ERG+ nuclei that PASSED the 75% rule
// in red (appearing magenta over DAPI) and only FINAL accepted PLVAP
// pixels in green. These are the objects used in the reported quantities.
selectWindow(ergPosMask);
run("Select None");
selectWindow(acceptedMask);
run("Select None");
selectWindow(dapiMask);
run("Select None");
run("Merge Channels...", "c1=[" + ergPosMask + "] c2=[" + acceptedMask + "] c3=[" + dapiMask + "] create keep");
rename(tripleAcceptedQC);
run("Select None");

print("3-color QC (blue=DAPI, red=thresholded ERG inside DAPI, green=ALL filtered PLVAP candidates): " + tripleCandidateQC);
print("3-color QC (blue=DAPI, magenta=ERG-positive nuclear ROIs passing 75%, green=ACCEPTED PLVAP): " + tripleAcceptedQC);
print("The three-color QC windows are for visual inspection only; quantitative results are unchanged.");

// ==================================================
// 9. Summary
// ==================================================

// ==================================================
// 8D. PLVAP-POSITIVE ENDOTHELIAL CELL COUNT (one count per ERG+ nucleus)
// ==================================================
// This applies the validated CAR4 cell-count definition to PLVAP,
// but retains the new DAPI+/ERG+ >=75% gate and the true calibrated 5-um radius.
// A PLVAP-positive "cell" is operationally an ERG+ nucleus whose NON-NUCLEAR
// 5-um local region has >=5% pixels with PLVAP intensity > a distant-background
// P99 threshold. This does NOT require >=5% of the nuclear area to be PLVAP+.
// Biological caveat: spatial association alone cannot assign all membrane PLVAP
// signal uniquely to one endothelial cell in confluent alveolar capillary networks.
// Unlike the object-count mask, this uses the P99 threshold, not mean+3SD.
// Both output types are kept separate and clearly labeled.

nonNuclearMask = prefix + "_NONNUCLEAR_valid_mask";
plvapCellBgExclusion = prefix + "_PLVAPcell_background_10um_EXCLUSION";
plvapCellBgValid = prefix + "_PLVAPcell_distant_BACKGROUND_valid";
plvapCellThresholdMask = prefix + "_PLVAPcell_P99_threshold_mask";
plvapCellCountNonNuc = prefix + "_PLVAPcell_P99_NONNUCLEAR_mask";
plvapCellPositiveMask = prefix + "_PLVAP_POSITIVE_ERG_NUCLEI_5pct_P99";
plvapCellQC = prefix + "_QC_DAPIblue_ERGred_PLVAPP99green";
plvapCellPerNucleusMask = prefix + "_TEMP_CELL_5um_zone";

// DAPI-positive regions must be excluded from the PLVAP cell-local fraction.
selectWindow(dapiMask);
run("Select None");
run("Duplicate...", "title=" + nonNuclearMask);
selectWindow(nonNuclearMask);
run("Invert");
run("Select None");

// Generate a true Euclidean 10-um exclusion region around ERG-positive nuclei
// to define distant PLVAP background, with no ROI Enlarge or Fill operations.
newImage(plvapCellBgExclusion, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
bgExclusionPixels = -floor(-(plvapCellBgExclusion_um / pixelSize_um)); // ImageJ macro compatible ceiling

if (nERGpos > 0) {
    roiManager("Reset");
    roiManager("Open", ergRoiZip);
    if (roiManager("count") != nERGpos)
        exit("CELL-COUNT QC ERROR: ERG-positive nuclear ROIs were not preserved.");
    for (ci=0; ci<nERGpos; ci++) {
        selectWindow(ergPosMask);
        roiManager("Select", ci);
        if (selectionType() == -1)
            exit("CELL-COUNT QC ERROR: invalid ERG nuclear ROI, index " + ci);
        Roi.getContainedPoints(nucX, nucY);
        selectWindow(plvapCellBgExclusion);
        for (cj=0; cj<nucX.length; cj++) {
            xx = nucX[cj]; yy = nucY[cj];
            for (ddy=-bgExclusionPixels; ddy<=bgExclusionPixels; ddy++) {
                qy = yy + ddy;
                if (qy>=0 && qy<imgHeight) {
                    for (ddx=-bgExclusionPixels; ddx<=bgExclusionPixels; ddx++) {
                        qx = xx + ddx;
                        r2 = (ddx*pixelWidth_um)*(ddx*pixelWidth_um) + (ddy*pixelHeight_um)*(ddy*pixelHeight_um);
                        if (qx>=0 && qx<imgWidth && r2<=plvapCellBgExclusion_um*plvapCellBgExclusion_um)
                            setPixel(qx,qy,255);
                    }
                }
            }
        }
    }
}
selectWindow(plvapCellBgExclusion);
run("Select None");
run("Invert"); // White = outside 10-um exclusion zone
selectWindow(nonNuclearMask);
run("Select None");
imageCalculator("AND create", nonNuclearMask, plvapCellBgExclusion);
rename(plvapCellBgValid);
run("Select None");

// Measure the histogram ONLY in the valid distant-background pixels.
// Crucially, this is a full 16-bit histogram (65,536 bins), not the legacy
// 256-bin, 8-bit histogram. With 16-bit data, bin index = intensity value.
selectWindow(plvapCellBgValid);
run("Select None");
getHistogram(validBgVals, validBgCounts, 256);
expectedBgPixels = validBgCounts[255];
run("Create Selection");
if (selectionType() == -1)
    exit("CELL-COUNT QC ERROR: no distant non-nuclear PLVAP background available.");
roiManager("Reset");
roiManager("Add");
selectWindow(plvapBgsub);
roiManager("Select", 0);
if (imageBits == 16) {
    getHistogram(histValues, histCounts, 65536);
    nHistBins = 65536;
} else {
    getHistogram(histValues, histCounts, 256);
    nHistBins = 256;
}
roiManager("Reset");
selectWindow(plvapBgsub);
run("Select None");
bgTotal = 0;
for (hi=0; hi<nHistBins; hi++) bgTotal += histCounts[hi];
if (bgTotal < 100)
    exit("CELL-COUNT QC ERROR: insufficient distant background pixels (<100).");
if (abs(bgTotal - expectedBgPixels) > 1)
    exit("CELL-COUNT QC ERROR: background selection included " + bgTotal + " pixels, but valid-background mask has " + expectedBgPixels + ". Check mask ROI fidelity.");
bgCum = 0;
plvapCellBgP99 = 0;
for (hi=0; hi<nHistBins; hi++) {
    bgCum += histCounts[hi];
    if (bgCum >= 0.99 * bgTotal) {
        plvapCellBgP99 = hi;
        break;
    }
}
print("PLVAP cell-count distant non-nuclear background: " + bgTotal + " pixels");
print("PLVAP cell-count distant-background P99 threshold (bg-subtracted intensity): " + plvapCellBgP99);

// Binary PLVAP-positive pixels for the cell-count rule.
// Distinct from the mean+3SD particle threshold used for PLVAP area/object counts.
selectWindow(plvapBgsub);
run("Select None");
run("Duplicate...", "title=" + plvapCellThresholdMask);
selectWindow(plvapCellThresholdMask);
getStatistics(plvapA, plvapMean, plvapMin, plvapMax);
if (plvapCellBgP99 >= plvapMax) {
    run("Select All"); run("Clear"); run("Select None");
    run("8-bit");
    print("WARNING: no PLVAP signal is above distant-background P99.");
} else {
    setThreshold(plvapCellBgP99 + 1, plvapMax);
    run("Convert to Mask");
    getStatistics(maskA, maskMean);
    if (maskMean > 127) run("Invert");
}
selectWindow(plvapCellThresholdMask); run("Select None");
selectWindow(nonNuclearMask); run("Select None");
imageCalculator("AND create", plvapCellThresholdMask, nonNuclearMask);
rename(plvapCellCountNonNuc);
run("Select None");

// Build a TRUE 5-um neighborhood for EACH nucleus. Use a temporary binary mask
// to prevent double-counting pixels reached from multiple pixels of one nucleus.
// Count local binary signal only outside the DAPI mask.
newImage(plvapCellPositiveMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
newImage(plvapCellPerNucleusMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
roiManager("Reset");
if (nERGpos > 0) roiManager("Open", ergRoiZip);
if (roiManager("count") != nERGpos)
    exit("CELL-COUNT QC ERROR: expected " + nERGpos + " ERG nuclei, but restored " + roiManager("count") + " ROIs.");

nPLVAPpositiveCells = 0;
nPLVAPcellNoValidLocalArea = 0;
for (ci=0; ci<nERGpos; ci++) {
    selectWindow(ergPosMask);
    roiManager("Select", ci);
    if (selectionType() == -1)
        exit("CELL-COUNT QC ERROR: invalid ERG nucleus index " + ci);
    Roi.getContainedPoints(nucX, nucY);

    selectWindow(plvapCellPerNucleusMask);
    run("Select All"); run("Clear"); run("Select None");
    for (cj=0; cj<nucX.length; cj++) {
        xx = nucX[cj]; yy = nucY[cj];
        for (ddy=-ergDilatePixels; ddy<=ergDilatePixels; ddy++) {
            qy = yy + ddy;
            if (qy>=0 && qy<imgHeight) {
                for (ddx=-ergDilatePixels; ddx<=ergDilatePixels; ddx++) {
                    qx = xx + ddx;
                    r2 = (ddx*pixelWidth_um)*(ddx*pixelWidth_um) + (ddy*pixelHeight_um)*(ddy*pixelHeight_um);
                    if (qx>=0 && qx<imgWidth && r2<=nearDistance_um*nearDistance_um)
                        setPixel(qx,qy,255);
                }
            }
        }
    }
    run("Select None");
    run("Create Selection");
    if (selectionType() == -1)
        exit("CELL-COUNT QC ERROR: no local zone for nucleus index " + ci);
    roiManager("Add");
    localROIindex = roiManager("count") - 1;
    selectWindow(nonNuclearMask);
    roiManager("Select", localROIindex);
    getHistogram(validVals, validCounts, 256);
    localValid = validCounts[255];
    selectWindow(plvapCellCountNonNuc);
    roiManager("Select", localROIindex);
    getHistogram(posVals, posCounts, 256);
    localPLVAP = posCounts[255];
    roiManager("Select", localROIindex);
    roiManager("Delete");

    localPLVAPFraction = 0;
    if (localValid > 0) localPLVAPFraction = localPLVAP / localValid;
    else nPLVAPcellNoValidLocalArea++;
    if (localValid > 0 && localPLVAPFraction >= plvapCellMinLocalFraction) {
        nPLVAPpositiveCells++;
        selectWindow(plvapCellPositiveMask);
        for (cj=0; cj<nucX.length; cj++)
            setPixel(nucX[cj], nucY[cj], 255);
    }
}
selectWindow(plvapCellPositiveMask); run("Select None");
// Keep plvapCellPerNucleusMask open for the DAPI+/PLVAP+ gate below.
if (nPLVAPcellNoValidLocalArea > 0)
    print("QC WARNING: " + nPLVAPcellNoValidLocalArea + " nuclei had no valid non-nuclear pixels in the local region.");
if (nPLVAPpositiveCells > nERGpos)
    exit("CELL-COUNT QC ERROR: PLVAP-positive cells exceed ERG-positive nuclei.");
plvapPositiveCellFraction = 0;
if (nERGpos > 0) plvapPositiveCellFraction = nPLVAPpositiveCells/nERGpos;

// Image for QC: blue DAPI, magenta DAPI+/ERG+ nuclei, green PLVAP above P99
// in non-nuclear pixels. The separate PLVAP-positive-nuclei mask shows which
// nuclei passed the >=5% gate.
selectWindow(ergPosMask); run("Select None");
selectWindow(plvapCellCountNonNuc); run("Select None");
selectWindow(dapiMask); run("Select None");
run("Merge Channels...", "c1=[" + ergPosMask + "] c2=[" + plvapCellCountNonNuc + "] c3=[" + dapiMask + "] create keep");
rename(plvapCellQC);
run("Select None");

// ==================================================
// 8E. DAPI+/PLVAP+ nuclei irrespective of ERG status
// ==================================================
// Apply EXACTLY the same PLVAP P99 threshold, 5-um local circular region,
// non-nuclear valid-pixel denominator, and >=5% rule as in section 8D.
// Nucleus is counted once, independently of the number of PLVAP fragments.
// The comparison is paired: every ERG+/PLVAP+ nucleus must also appear here.
// This is an operational spatial-association count, not definite cell membrane
// assignment or proof that a PLVAP-positive cell is ERG negative.

dapiPLVAPposMask = prefix + "_DAPI_PLVAP_POSITIVE_NUCLEI_5pct_P99";
dapiPLVAPQC = prefix + "_QC_DAPIblue_PLVAPposNUCLEIred_PLVAPP99green";
newImage(dapiPLVAPposMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");

roiManager("Reset");
roiManager("Open", dapiRoiZip);
if (roiManager("count") != nDAPI)
    exit("DAPI/PLVAP CELL-COUNT QC ERROR: expected " + nDAPI + " DAPI nuclear ROIs, restored " + roiManager("count") + ".");

nDAPIPLVAPpositiveCells = 0;
nDAPIPLVAPnoValidLocalArea = 0;
for (ci=0; ci<nDAPI; ci++) {
    selectWindow(dapiMask);
    roiManager("Select", ci);
    if (selectionType() == -1)
        exit("DAPI/PLVAP QC ERROR: invalid DAPI nucleus at ROI index " + ci);
    Roi.getContainedPoints(nucX, nucY);
    if (nucX.length == 0)
        exit("DAPI/PLVAP QC ERROR: empty DAPI nucleus at ROI index " + ci);

    selectWindow(plvapCellPerNucleusMask);
    run("Select All"); run("Clear"); run("Select None");
    for (cj=0; cj<nucX.length; cj++) {
        xx = nucX[cj]; yy = nucY[cj];
        for (ddy=-ergDilatePixels; ddy<=ergDilatePixels; ddy++) {
            qy = yy + ddy;
            if (qy>=0 && qy<imgHeight) {
                for (ddx=-ergDilatePixels; ddx<=ergDilatePixels; ddx++) {
                    qx = xx + ddx;
                    r2 = (ddx*pixelWidth_um)*(ddx*pixelWidth_um) + (ddy*pixelHeight_um)*(ddy*pixelHeight_um);
                    if (qx>=0 && qx<imgWidth && r2<=nearDistance_um*nearDistance_um)
                        setPixel(qx,qy,255);
                }
            }
        }
    }
    run("Select None");
    run("Create Selection");
    if (selectionType() == -1)
        exit("DAPI/PLVAP QC ERROR: no local zone for DAPI ROI index " + ci);
    roiManager("Add");
    localROIindex = roiManager("count") - 1;
    selectWindow(nonNuclearMask);
    roiManager("Select", localROIindex);
    getHistogram(validVals, validCounts, 256);
    localValid = validCounts[255];
    selectWindow(plvapCellCountNonNuc);
    roiManager("Select", localROIindex);
    getHistogram(posVals, posCounts, 256);
    localPLVAP = posCounts[255];
    roiManager("Select", localROIindex);
    roiManager("Delete");

    localPLVAPFraction = 0;
    if (localValid > 0) localPLVAPFraction = localPLVAP / localValid;
    else nDAPIPLVAPnoValidLocalArea++;
    if (localValid > 0 && localPLVAPFraction >= plvapCellMinLocalFraction) {
        nDAPIPLVAPpositiveCells++;
        selectWindow(dapiPLVAPposMask);
        for (cj=0; cj<nucX.length; cj++)
            setPixel(nucX[cj], nucY[cj], 255);
    }
}
selectWindow(dapiPLVAPposMask); run("Select None");
selectWindow(plvapCellPerNucleusMask); close();

if (nDAPIPLVAPnoValidLocalArea > 0)
    print("QC WARNING: " + nDAPIPLVAPnoValidLocalArea + " DAPI nuclei lacked a non-nuclear local region for PLVAP testing.");
if (nDAPIPLVAPpositiveCells > nDAPI)
    exit("DAPI/PLVAP CELL-COUNT QC ERROR: positive count exceeds number of DAPI nuclei.");
if (nPLVAPpositiveCells > nDAPIPLVAPpositiveCells)
    exit("DAPI/PLVAP CELL-COUNT QC ERROR: ERG+/PLVAP+ cell count exceeds total DAPI+/PLVAP+ count, despite identical rules. Check ROI and mask alignment.");

dapiPLVAPfraction = nDAPIPLVAPpositiveCells / nDAPI;
nDAPIPLVAP_ERGnegative = nDAPIPLVAPpositiveCells - nPLVAPpositiveCells;

selectWindow(dapiPLVAPposMask); run("Select None");
selectWindow(plvapCellCountNonNuc); run("Select None");
selectWindow(dapiMask); run("Select None");
run("Merge Channels...", "c1=[" + dapiPLVAPposMask + "] c2=[" + plvapCellCountNonNuc + "] c3=[" + dapiMask + "] create keep");
rename(dapiPLVAPQC);
run("Select None");

File.delete(dapiRoiZip);
File.delete(ergRoiZip);
print("DAPI+/PLVAP+ NUCLEI, regardless of ERG (5% local/P99): " + nDAPIPLVAPpositiveCells + " / " + nDAPI + " = " + dapiPLVAPfraction);
print("DAPI+/PLVAP+ nuclei that FAILED the ERG gate: " + nDAPIPLVAP_ERGnegative + " (not proven ERG-negative cells).");
print("New DAPI/PLVAP cell-count QC images: " + dapiPLVAPposMask + " and " + dapiPLVAPQC);
print("PLVAP-POSITIVE ERG+ ENDOTHELIAL CELLS (cell-level 5%/P99 rule): " + nPLVAPpositiveCells + " / " + nERGpos + " = " + plvapPositiveCellFraction);
print("Cell rule: >= " + 100*plvapCellMinLocalFraction + "% of non-nuclear pixels within circular " + nearDistance_um + " um region > same-image distant-background P99.");
print("These are cell-level PLVAP-associated estimates, NOT counts of PLVAP mask fragments.");
print("Cell count QC images: " + plvapCellPositiveMask + " and " + plvapCellQC);


// ==================================================
// 9. Summary
// ==================================================

print("------ SUMMARY ------");
print("DAPI nuclei detected: " + nDAPI);
print("ERG-positive nuclei detected: " + nERGpos);
print("ERG-positive fraction of DAPI nuclei: " + ergPosFraction);
print("ERG threshold method: " + ergThresholdMethod);
print("ERG SD multiplier: " + ergSDmultiplier);
print("ERG nuclear fraction threshold: " + ergNuclearFractionThreshold);
print("ERG expansion distance, um: " + nearDistance_um);
print("ERG_near_zone positive fraction of image: " + ergZoneFraction);
print("PLVAP-associated ERG-positive nuclei: " + nPLVAPassocERG);
print("PLVAP-associated fraction of ERG-positive nuclei: " + plvapAssocERGFraction);
print("DAPI+/PLVAP+ nuclei, regardless of ERG (5% local / P99): " + nDAPIPLVAPpositiveCells);
print("DAPI+/PLVAP+ fraction of DAPI nuclei (5% local / P99): " + dapiPLVAPfraction);
print("DAPI+/PLVAP+ nuclei without >=75% ERG overlap: " + nDAPIPLVAP_ERGnegative);
print("DAPI+/ERG+/PLVAP+ nuclei (5% local / P99): " + nPLVAPpositiveCells);
print("PLVAP+ fraction of ERG-positive nuclei (5% local / P99): " + plvapPositiveCellFraction);
print("PLVAP cell rule: >=5% local non-nuclear pixels within 5 um above distant-background P99=" + plvapCellBgP99);
print("PLVAP threshold method: " + plvapThresholdMethod);
print("PLVAP candidate objects before ERG-distance filtering: " + nPLVAPcandidate);
print("ERG-distance-associated PLVAP objects: " + nAccepted);
print("Total ERG-distance-associated PLVAP area, um2: " + totalArea);
print("Total raw PLVAP integrated density: " + totalRawIntegratedDensity);
print("Total background-subtracted PLVAP integrated density: " + totalBgSubIntegratedDensity);

if (nERGpos > 0) {
    print("ERG-distance-associated PLVAP objects per ERG-positive nucleus: " + nAccepted / nERGpos);
    print("Background-subtracted PLVAP integrated density per ERG-positive nucleus: " + totalBgSubIntegratedDensity / nERGpos);
}

print("---------------------");
// ==================================================
// 10. THREE FINAL COLORED QC IMAGES (DISPLAY ONLY)
// No changes to any counts, gating, thresholds or distance tests.
// ImageJ Merge Channels: c1=red, c2=green, c3=blue.
// Because selected nuclear pixels overlap blue DAPI, they display magenta.
// ==================================================
// Short per-run identifier prevents accidental reuse of an existing QC window.
qcRunTag = d2s(getTime()%1000000, 0);
qcDAPIName = "QC1_PLVAP_DAPI_BLUE_n" + nDAPI + "_" + qcRunTag;
qcERGName = "QC2_PLVAP_DAPIblue_ERGmagenta_n" + nERGpos + "_" + qcRunTag;
qcTripleName = "QC3_PLVAP_TRIPLEmagenta_PLVAPgreen_n" + nPLVAPpositiveCells + "_" + qcRunTag;

// Verify segmentation masks and above-P99 marker map are available.
qcTitlesBefore = getList("image.titles");
qcDAPIok = 0; qcERGok = 0; qcTripleOK = 0; qcMarkerOK = 0;
for (qi=0; qi<qcTitlesBefore.length; qi++) {
    if (qcTitlesBefore[qi] == dapiDetectedQCMask) qcDAPIok = 1;
    if (qcTitlesBefore[qi] == ergPosMask) qcERGok = 1;
    if (qcTitlesBefore[qi] == plvapCellPositiveMask) qcTripleOK = 1;
    if (qcTitlesBefore[qi] == plvapCellCountNonNuc) qcMarkerOK = 1;
}
if (qcDAPIok == 0 || qcERGok == 0 || qcTripleOK == 0 || qcMarkerOK == 0)
    exit("FINAL COLOR QC ERROR: source QC masks missing. Analysis log retained; no intermediate files closed.");

// QC1: Only actually counted DAPI nuclei, blue on black. Never alter source mask.
selectWindow(dapiDetectedQCMask); run("Select None");
run("Duplicate...", "title=" + qcDAPIName);
selectWindow(qcDAPIName);
run("Select None");
run("Blue");

// QC2: All DAPI nuclei in blue; nuclei that passed 75% DAPI/ERG in magenta.
selectWindow(ergPosMask); run("Select None");
selectWindow(dapiDetectedQCMask); run("Select None");
run("Merge Channels...", "c1=[" + ergPosMask + "] c3=[" + dapiDetectedQCMask + "] create keep");
rename(qcERGName);
run("Select None");

// QC3: Counted triple-positive nuclei are magenta, DAPI blue; bright
// PLVAP pixels above background P99 are green as CONTEXT, not cell counts.
// Only magenta nuclear regions in QC3 count towards the triple-positive total.
selectWindow(plvapCellPositiveMask); run("Select None");
selectWindow(plvapCellCountNonNuc); run("Select None");
selectWindow(dapiDetectedQCMask); run("Select None");
run("Merge Channels...", "c1=[" + plvapCellPositiveMask + "] c2=[" + plvapCellCountNonNuc + "] c3=[" + dapiDetectedQCMask + "] create keep");
rename(qcTripleName);
run("Select None");

// Keep only the 3 colored QC windows generated this run and the input TIFF.
// Close only images carrying this run's unique prefix; other runs are untouched.
qcList = getList("image.titles");
closedQC = 0;
for (qi=0; qi<qcList.length; qi++) {
    if (indexOf(qcList[qi], prefix + "_") == 0) {
        selectWindow(qcList[qi]);
        run("Select None");
        close();
        closedQC++;
    }
}
print("---------------------");
print("THREE FINAL COLORED QC WINDOWS:");
print("QC1: Blue counted DAPI nuclei n=" + nDAPI + " (" + qcDAPIName + ")");
print("QC2: Blue DAPI + magenta DAPI/ERG nuclei n=" + nERGpos + " (" + qcERGName + ")");
print("QC3: Blue DAPI + magenta triple-positive nuclei n=" + nPLVAPpositiveCells + " + green above-P99 PLVAP context (" + qcTripleName + ")");
print("Closed " + closedQC + " intermediate images. Original TIFF remains open.");
print("GREEN represents PLVAP fluorescence above P99, NOT separate cells; only MAGENTA nuclei in QC3 are counted as triple-positive.");
