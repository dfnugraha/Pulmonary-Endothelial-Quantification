// CAR4 cell counting: final QC shows ONLY DAPI, DAPI+/ERG+, DAPI+/ERG+/CAR4+ nuclei
// Adapted for OlyVIA 3-plane 8/16-bit raw TIFF, including 1539 x 1539 images
// Uses the validated 75% nuclear-overlap / 5-um proximity workflow.
// IMPORTANT: This is a proximity-based count. The older CAR4 macro used a
// distinct >=5% non-nuclear-local-area / background P99 definition.
// Those definitions must not be pooled or silently compared.
// Reads physical pixel size from TIFF calibration. Does not assume a 512 px field.
//
// Logic:
// DAPI nuclei -> ERG-positive nuclei -> exact circular 5-um proximity zone -> CAR4 candidate mask
// -> accepted CAR4 = CAR4_mask AND ERG_near_zone. QC includes all gates.
//
// Expected raw Olympus VSI channel order for the supplied CAR4 acquisitions:
// C1 = DAPI
// C2 = CAR4 (verify actual fluorophore in the staining protocol)
// C3 = ERG (verify actual fluorophore in the staining protocol)
//
// IMPORTANT: display pseudocolors in OlyVIA do not change this raw channel order.
//
// Main fixes:
// Final association now checks accepted CAR4 pixels within a fixed
// neighborhood of each original DAPI+/ERG+ nuclear ROI, without using
// ROI Enlarge (which failed for ROI index 2).
// The previous version left the last nucleus ROI selected.
// In Fiji, Duplicate can duplicate only the active selection, not the full image.
// This made ERG_near_zone wrong. This version clears selections before duplicating masks.

requires("1.53");

orig = getTitle();
prefix = "ERG_CAR4_CellCount_75pct_" + d2s(getTime(), 0);

// ================= USER SETTINGS =================

cERG   = 3;
cDAPI  = 1;
cCAR4 = 2;

// Pixel size will be read from ImageJ calibration in the TIFF.

// Exact Euclidean 5-um distance from an accepted DAPI+/ERG+ nuclear pixel
nearDistance_um = 5;

// Background subtraction
// Spatially scaled from full-res 30/50/100 px at 0.325012 um/px
rollingBallDAPI  = 10;
rollingBallERG   = 17;
rollingBallCAR4 = 33;

// Spatially scaled Gaussian sigmas (full-res: 1/1/0.3 px)
gaussianSigmaDAPI  = 0.33;
gaussianSigmaERG   = 0.33;
gaussianSigmaCAR4 = 0.10;

// ERG threshold carried forward from the supplied 512 x 512 macro.
// The 488/ERG channel contains substantial diffuse septal signal.
// Strict test: mean + 2.5 SD after background subtraction, and a DAPI nucleus
// must have >=75% overlap with above-threshold ERG pixels. This is not a
// biological truth: inspect nuclear QC masks, and compare with previous settings.
// Both criteria are configurable below.
dapiThresholdMethod  = "Otsu";
ergThresholdMethod   = "Manual_SD";
ergSDmultiplier         = 2.5;      // stricter ERG threshold, previously 2 SD
car4ThresholdMethod = "Manual_SD";  // stricter: mean + SD multiplier threshold

// DAPI nucleus filtering, um^2
dapiMinArea = 15;
dapiMaxArea = 300;

// Strict ERG-positive nucleus definition (same as the validated PLVAP workflow)
ergNuclearFractionThreshold = 0.75;  // >=75% of DAPI nucleus must exceed ERG intensity threshold

// CAR4 object filtering, um^2
// Unlike punctate PLVAP, CAR4 may form extended alveolar-septal structures.
// Preserve the older CAR4 macro convention: NO maximum object-size cutoff.
car4MinArea = 0.20;
car4MaxArea = "Infinity";

// CAR4-specific intensity threshold carried over from the earlier CAR4 macro.
// Note: CAR4 morphology differs from punctate PLVAP.
car4SDmultiplier = 3.0;

// Original CAR4final(2).ijm per-cell classification rule, adapted to 16-bit TIFF:
// >=5% of non-nuclear pixels within 5 um of each ERG+ nucleus must be brighter
// than the 99th percentile of same-image distant non-nuclear CAR4 background.
car4CellMinLocalFraction = 0.05;
car4CellBgExclusion_um = 10;

// Preserve extended CAR4-positive structures; do not filter by circularity.
car4CircularityMin = 0.00;

dapiWatershed = true;
// Keep Despeckle off to preserve CAR4 signal morphology.
car4DoDespeckle = false;

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
rollingBallCAR4 = round(rollingBallCAR4 * referencePixelSize_um / pixelSize_um);
if (rollingBallCAR4 < 1) rollingBallCAR4 = 1;
gaussianSigmaDAPI = gaussianSigmaDAPI * referencePixelSize_um / pixelSize_um;
gaussianSigmaERG = gaussianSigmaERG * referencePixelSize_um / pixelSize_um;
gaussianSigmaCAR4 = gaussianSigmaCAR4 * referencePixelSize_um / pixelSize_um;

run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
run("Set Measurements...", "area mean min max integrated redirect=None decimal=3");
setOption("BlackBackground", true);
setBackgroundColor(0, 0, 0);
setForegroundColor(255, 255, 255);

ergDilatePixels = -floor(-(nearDistance_um / pixelSize_um));  // ImageJ macro compatible ceiling; actual membership uses physical Euclidean distance

print("\\Clear");
print("Image: " + orig);
print("Unique run prefix: " + prefix);
print("Image size: " + imgWidth + " x " + imgHeight);
print("Channels: " + imgChannels + ", Slices: " + imgSlices + ", Frames: " + imgFrames);
print("Input: " + imgWidth + " x " + imgHeight + ", C=" + imgChannels + ", " + imageBits + "-bit");
print("Calibrated field: " + fieldWidth_um + " x " + fieldHeight_um + " um");
print("Pixel size: " + pixelSize_um + " um/pixel");
print("Raw channel order: C1 DAPI / C2 CAR4 / C3 ERG (VERIFY!)");
print("ERG expansion distance: " + nearDistance_um + " um");
print("ERG distance-test search radius in pixels: " + ergDilatePixels);
print("ERG proximity uses physical Euclidean distance <= " + nearDistance_um + " um (not a square Dilate neighborhood).");
print("ERG association uses direct pixel-neighborhood tests; no ROI enlargement.");

// Window names
ergRaw = prefix + "_ERG_raw";
dapiRaw = prefix + "_DAPI_raw";
car4Raw = prefix + "_CAR4_raw";
dapiMask = prefix + "_DAPI_mask";
ergThresh = prefix + "_ERG_thresh_mask";
ergDapiOverlapMask = prefix + "_ERG_above_threshold_INSIDE_DAPI_QC";
ergTwoColorQC = prefix + "_DAPI_blue_ERG_green_raw_QC";
ergNuclearTwoColorQC = prefix + "_DAPI_blue_ERG_nuclear_red_QC";
ergPosMask = prefix + "_ERGpos_nuclei_mask";
ergZone = prefix + "_ERG_near_zone_TRUE_circular_5um";
car4Bgsub = prefix + "_CAR4_bgsub";
car4Mask = prefix + "_CAR4_threshold_mask";
car4FilteredMask = prefix + "_CAR4_candidate_filtered_mask";
acceptedPreMask = prefix + "_Accepted_CAR4_pre_filter";
acceptedMask = prefix + "_Accepted_CAR4_mask";
acceptedRawMap = prefix + "_Accepted_CAR4_raw_intensity_map";
acceptedBgsubMap = prefix + "_Accepted_CAR4_bgsub_intensity_map";
car4AssocERGMask = prefix + "_CAR4_associated_ERGpos_nuclei_mask";
car4NearForERG = prefix + "_CAR4_near_for_ERG_count";
associationOverlapQC = prefix + "_CAR4near_AND_ERGpos_QC";
tripleCandidateQC = prefix + "_DAPI_blue_ERG_nuclear_red_CAR4_candidates_green_QC";
tripleAcceptedQC = prefix + "_DAPI_blue_ERGpos_red_ACCEPTED_CAR4_green_QC";
tripleGateQC = prefix + "_QC_accepted_ERGmagenta_ALL_CAR4_green_DAPIblue";
tripleDistanceQC = prefix + "_QC_5umZONEred_CAR4_candidatesGREEN_accepted_ERGblue";

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
Stack.setChannel(cCAR4);
run("Duplicate...", "title=" + car4Raw + " duplicate channels=" + cCAR4);

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

// Preserve ALL segmented DAPI nuclei for the new DAPI+/CAR4+ count.
// The ERG-positive ROI set is saved separately later.
dapiRoiZip = getDirectory("temp") + prefix + "_ALL_DAPI_nuclear_ROIs.zip";
roiManager("Save", dapiRoiZip);

// QC DISPLAY ONLY: Build a mask from the nuclear ROIs actually counted.
// The DAPI threshold image can contain rejected spots; it is not a count mask.
dapiDetectedQCMask = prefix + "_DAPI_COUNTED_NUCLEI_mask";
newImage(dapiDetectedQCMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
for (di=0; di<nDAPI; di++) {
    selectWindow(dapiMask);
    roiManager("Select", di);
    if (selectionType() == -1) exit("DAPI COLOR QC ERROR: missing counted DAPI ROI " + di);
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
// The ERG threshold is applied to the actual C3 ERG channel.
// Retain only pixels above mean + 2.5 SD of the background-subtracted ERG channel.
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

// ERG double-positive threshold is fixed at 75%, with QC of mask morphology.
// The macro now continues quantification and reports a warning only.
if (ergPosFraction > 0.80) {
    print("WARNING: More than 80% of DAPI nuclei are classified as ERG-positive.");
    print("Quantification will continue. Inspect ERG_thresh_mask and ERGpos_nuclei_mask for QC.");
}

// Preserve exact ERG-positive DAPI ROIs before the ROI Manager is reused
// for CAR4 particle filtering. This lets us later count ERG+ nuclei that
// have accepted CAR4 signal within the same 5-um association distance.
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

// Build a TRUE circular 5-um zone from accepted ERG-positive nuclei.
// Repeated binary Dilate produces a square/Chebyshev 5-pixel neighborhood
// and allows diagonally separated pixels beyond 5 um. Instead use calibrated
// x/y distances and save each valid pixel directly. No ROI Enlarge or Fill.
newImage(ergZone, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");

if (nERGpos > 0) {
    roiManager("Reset");
    roiManager("Open", ergRoiZip);
    if (roiManager("count") != nERGpos)
        exit("ROI QC FAILURE: the saved ERG-positive nuclei could not be reloaded for 5-um distance zone.");

    for (i=0; i<nERGpos; i++) {
        selectWindow(ergPosMask);
        roiManager("Select", i);
        if (selectionType() == -1)
            exit("ROI QC FAILURE: missing ERG-positive nucleus ROI in distance-zone generation, index " + i);
        Roi.getContainedPoints(nucleusX, nucleusY);
        if (nucleusX.length == 0)
            exit("ROI QC FAILURE: empty ERG-positive nucleus ROI in distance-zone generation, index " + i);

        selectWindow(ergZone);
        for (j=0; j<nucleusX.length; j++) {
            nx = nucleusX[j];
            ny = nucleusY[j];
            for (dy=-ergDilatePixels; dy<=ergDilatePixels; dy++) {
                py = ny + dy;
                if (py>=0 && py<imgHeight) {
                    for (dx=-ergDilatePixels; dx<=ergDilatePixels; dx++) {
                        px = nx + dx;
                        distanceSquared_um2 = (dx * pixelWidth_um) * (dx * pixelWidth_um) +
                            (dy * pixelHeight_um) * (dy * pixelHeight_um);
                        if (px>=0 && px<imgWidth && distanceSquared_um2 <= nearDistance_um * nearDistance_um)
                            setPixel(px, py, 255);
                    }
                }
            }
        }
    }
}
selectWindow(ergZone);
run("Select None");
getStatistics(areaZoneCheck, meanZoneCheck);
ergZoneFraction = meanZoneCheck / 255;
print("ERG_near_zone positive fraction of image (Euclidean distance): " + ergZoneFraction);

if (ergZoneFraction > 0.80) {
    print("WARNING: ERG_near_zone covers more than 80% of the image.");
    print("Reduce ERG-positive nuclei count or reduce nearDistance_um.");
}

// ==================================================
// 6. CAR4 mask
// ==================================================

selectWindow(car4Raw);
run("Select None");
run("Duplicate...", "title=" + car4Bgsub);

selectWindow(car4Bgsub);
run("Subtract Background...", "rolling=" + rollingBallCAR4);

selectWindow(car4Bgsub);
run("Select None");
run("Duplicate...", "title=" + car4Mask);

selectWindow(car4Mask);

// Preserve CAR4 signal morphology without strong additional smoothing.
run("Gaussian Blur...", "sigma=" + gaussianSigmaCAR4);

// CAR4 candidate threshold: background-subtracted mean + 3 SD,
// with no upper area limit or circularity exclusion.
getStatistics(areaCAR4stat, meanCAR4stat, minCAR4stat, maxCAR4stat, stdCAR4stat);
car4LowerThreshold = meanCAR4stat + car4SDmultiplier * stdCAR4stat;
setThreshold(car4LowerThreshold, maxCAR4stat);
print("CAR4 manual threshold lower value: " + car4LowerThreshold);
print("CAR4 bgsub mean: " + meanCAR4stat + ", SD: " + stdCAR4stat + ", max: " + maxCAR4stat);
run("Convert to Mask");

// Make sure CAR4 signal is white and background is black.
getStatistics(areaCAR4check, meanCAR4check);
if (meanCAR4check > 127) {
    run("Invert");
}

if (car4DoDespeckle == true) {
    run("Despeckle");
}

run("Select None");
roiManager("Reset");
// Ask ImageJ to make the mask directly from particles passing the filters.
// Single-pixel particles can have ROI types for which run("Fill") fails.
// The Analyze Particles mask includes these particles without ROI filling.
run("Analyze Particles...",
    "size=" + car4MinArea + "-" + car4MaxArea +
    " circularity=" + car4CircularityMin + "-1.00 show=Masks exclude add");

nCAR4candidate = roiManager("count");
if (getTitle() == car4Mask)
    exit("Could not create CAR4 particle mask from Analyze Particles. Verify Fiji supports show=Masks.");
rename(car4FilteredMask);
run("Select None");

selectWindow(car4FilteredMask);
run("Select None");
getStatistics(areaCAR4mask, meanCAR4mask);
car4MaskFraction = meanCAR4mask / 255;

print("Total CAR4 candidate objects before ERG-distance filtering: " + nCAR4candidate);
print("Filtered CAR4 candidate-mask positive fraction of image: " + car4MaskFraction);

// ==================================================
// 7. Accepted CAR4 = FILTERED CAR4 candidates AND ERG_near_zone
// ==================================================

selectWindow(car4FilteredMask);
run("Select None");
selectWindow(ergZone);
run("Select None");

imageCalculator("AND create", car4FilteredMask, ergZone);
rename(acceptedPreMask);

selectWindow(acceptedPreMask);
run("Select None");

getStatistics(areaAcceptedPre, meanAcceptedPre);
acceptedPreFraction = meanAcceptedPre / 255;
print("Accepted CAR4 mask positive fraction before final particle filtering: " + acceptedPreFraction);

if (car4DoDespeckle == true) {
    run("Despeckle");
}

run("Select None");
roiManager("Reset");
// Use the particle-mask output rather than ROI-by-ROI filling (which
// can fail for single-pixel accepted objects).
run("Analyze Particles...",
    "size=" + car4MinArea + "-" + car4MaxArea +
    " circularity=" + car4CircularityMin + "-1.00 show=Masks exclude add");

nAccepted = roiManager("count");
if (getTitle() == acceptedPreMask)
    exit("Could not create accepted CAR4 particle mask from Analyze Particles.");
rename(acceptedMask);
run("Select None");

selectWindow(acceptedMask);
run("Select None");
getStatistics(areaAcceptedMask, meanAcceptedMask);
acceptedMaskFraction = meanAcceptedMask / 255;

print("ERG-distance-associated CAR4 objects: " + nAccepted);
print("Final accepted CAR4 mask positive fraction: " + acceptedMaskFraction);

// ==================================================
// 7B. Create intensity-map QC images
// These show raw/bg-subtracted CAR4 intensity only inside accepted CAR4 objects.
// Everything outside Accepted_CAR4_mask is cleared to black.
// ==================================================

// If there are no accepted objects, a zero-mask has no area selection.
// Produce empty intensity maps rather than running Create Selection on it.
selectWindow(car4Raw);
run("Select None");
run("Duplicate...", "title=" + acceptedRawMap);
selectWindow(car4Bgsub);
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
    print("No accepted CAR4 objects: intensity QC maps set to zero.");
}

print("Created intensity QC maps:");
print("1. " + acceptedRawMap);
print("2. " + acceptedBgsubMap);

if (ergZoneFraction > 0.95 && car4MaskFraction > 0.001 && nAccepted == 0) {
    print("TECHNICAL WARNING: ERG_near_zone is almost full image, CAR4 mask is positive, but accepted CAR4 is 0.");
    print("This should not happen. Close old mask windows and rerun, or check mask polarity.");
}

// ==================================================
// 8. Measure accepted CAR4 intensity
// ==================================================

run("Clear Results");

totalArea = 0;
totalRawIntegratedDensity = 0;
totalBgSubIntegratedDensity = 0;

for (i = 0; i < nAccepted; i++) {
    selectWindow(car4Raw);
    roiManager("Select", i);
    getStatistics(areaRaw, meanRaw, minRaw, maxRaw, stdRaw);
    rawIntegratedDensity = areaRaw * meanRaw;

    selectWindow(car4Bgsub);
    roiManager("Select", i);
    getStatistics(areaBg, meanBg, minBg, maxBg, stdBg);
    bgSubIntegratedDensity = areaBg * meanBg;

    setResult("Accepted_CAR4_object_ID", i, i + 1);
    setResult("CAR4_area_um2", i, areaRaw);
    setResult("Mean_CAR4_raw", i, meanRaw);
    setResult("Integrated_density_raw", i, rawIntegratedDensity);
    setResult("Mean_CAR4_bgsub", i, meanBg);
    setResult("Integrated_density_bgsub", i, bgSubIntegratedDensity);

    totalArea += areaRaw;
    totalRawIntegratedDensity += rawIntegratedDensity;
    totalBgSubIntegratedDensity += bgSubIntegratedDensity;
}

updateResults();

// ==================================================
// 8B. CAR4-associated DAPI+/ERG+ nuclei -- NO ROI enlargement
// ==================================================
// A nucleus is associated when at least one final accepted CAR4-mask
// pixel is within the TRUE calibrated Euclidean 5-um distance from any
// original DAPI+/ERG+ nuclear pixel. Matches circular ERG_near_zone above.
// This corrects the former square neighborhood (which could include pixels
// >5 um away diagonally).
//
// No Enlarge, Restore Selection, or per-particle ROI dilation is used.
// The stored ROIs are used only to find nuclear pixels and to draw the
// positive nucleus into the final QC mask. If there is a mismatch,
// report a technical failure rather than treating it as zero biology.

nCAR4assocERG = 0;
car4AssocERGFraction = 0;
associationOverlapQC = prefix + "_Accepted_CAR4_AND_ERG_near_zone_QC";

newImage(car4AssocERGMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");

// QC: the accepted particle mask must actually overlap the nuclear near-zone.
selectWindow(acceptedMask);
run("Select None");
getStatistics(aAcceptedQC, mAcceptedQC, minAcceptedQC, maxAcceptedQC);
print("Final CAR4 mask pixel range: " + minAcceptedQC + " to " + maxAcceptedQC);
if (nAccepted > 0 && (maxAcceptedQC == 0 || minAcceptedQC != 0))
    exit("QC FAILURE: positive CAR4 objects were reported, but the final mask is missing/invalid.");

selectWindow(ergZone);
run("Select None");
selectWindow(acceptedMask);
run("Select None");
imageCalculator("AND create", acceptedMask, ergZone);
rename(associationOverlapQC);
run("Select None");
getStatistics(aAssocQC, mAssocQC, minAssocQC, maxAssocQC);
print("Accepted CAR4 overlapping ERG near-zone (mask fraction): " + (mAssocQC / 255));
if (nAccepted > 0 && maxAssocQC == 0)
    exit("QC FAILURE: accepted CAR4 particles do not overlap the ERG near-zone. Inspect mask polarity and image alignment.");

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
                        distanceSquared_um2 = (dx * pixelWidth_um) * (dx * pixelWidth_um) +
                            (dy * pixelHeight_um) * (dy * pixelHeight_um);
                        if (px>=0 && px<imgWidth && distanceSquared_um2 <= nearDistance_um * nearDistance_um && getPixel(px, py)>0)
                            associated = 1;
                    }
                }
            }
        }

        if (associated==1) {
            nCAR4assocERG++;
            // Draw from the pixel coordinates retrieved from the ORIGINAL nucleus ROI.
            // Do not transfer an ROI selection between windows: some Fiji ROI types
            // lose the selection when switching images, causing run("Fill") to fail.
            selectWindow(car4AssocERGMask);
            for (k=0; k<nucleusX.length; k++)
                setPixel(nucleusX[k], nucleusY[k], 255);
        }
    }
    selectWindow(car4AssocERGMask);
    run("Select None");
    car4AssocERGFraction = nCAR4assocERG / nERGpos;
    print("CAR4-proximity-associated DAPI+/ERG+ nuclei (>=1 accepted CAR4 pixel within 5 um): " + nCAR4assocERG + " / " + nERGpos + " = " + car4AssocERGFraction);
    print("Association searched original nuclear pixels within a calibrated circular " + nearDistance_um + " um radius, matching ERG_near_zone.");

    if (nAccepted > 0 && nCAR4assocERG==0)
        exit("ASSOCIATION QC FAILURE: CAR4 objects overlap the ERG-near-zone, but no original nucleus has a nearby CAR4 pixel. Inspect the nucleus mask/ROIs and binary masks; do not report zero.");
    if (nAccepted==0 && nCAR4assocERG>0)
        exit("ASSOCIATION QC FAILURE: detected associated nuclei with no accepted CAR4 objects.");

    // Preserve ROI zip for the separate CAR4-positive cell count below.
}

// ==================================================
// 8B-QC. Stage-by-stage proximity diagnostics, without changing counts
// ==================================================
// Pixel counts distinguish candidate CAR4 lost to the proximity gate
// versus CAR4 lost later to particle-area filtering. These are pixel
// counts, NOT counts of biological CAR4-positive cells.
totalCandidatePixels = round(car4MaskFraction * imgWidth * imgHeight);
insideDistancePixels = round(acceptedPreFraction * imgWidth * imgHeight);
finalAcceptedPixels = round(acceptedMaskFraction * imgWidth * imgHeight);
print("DISTANCE QC: CAR4 candidate-positive pixels = " + totalCandidatePixels);
print("DISTANCE QC: CAR4 pixels INSIDE true 5-um ERG zone, before final particle filter = " + insideDistancePixels);
print("DISTANCE QC: CAR4 pixels OUTSIDE true 5-um ERG zone = " + (totalCandidatePixels - insideDistancePixels));
print("DISTANCE QC: CAR4 pixels removed after proximity gate by final particle filter = " + (insideDistancePixels - finalAcceptedPixels));
print("A candidate CAR4 pixel outside the 5-um zone will remain green in the distance QC overlay, while one within the zone appears yellow (red + green).");

// ==================================================
// 8C. Three-color QC overlays (VISUALIZATION ONLY; new 5-um gate above affects counts)
// ==================================================
// These are display composites of the ALREADY CALCULATED masks.
// These overlays do not change counts. The proximity algorithm above was corrected to a true 5-um Euclidean radius.
// Fiji's Merge Channels mapping is: c1=RED, c2=GREEN, c3=BLUE.
// DAPI blue + ERG red overlap looks magenta on screen.
//
// Candidate QC: Shows threshold-positive ERG pixels within DAPI nuclei
// in red and all particle-filtered CAR4 candidates in green.
// Note that red pixels can occur in nuclei NOT passing the 75% ERG gate.
selectWindow(ergDapiOverlapMask);
run("Select None");
selectWindow(car4FilteredMask);
run("Select None");
selectWindow(dapiMask);
run("Select None");
run("Merge Channels...", "c1=[" + ergDapiOverlapMask + "] c2=[" + car4FilteredMask + "] c3=[" + dapiMask + "] create keep");
rename(tripleCandidateQC);
run("Select None");

// Accepted QC: Marks complete DAPI+/ERG+ nuclei that PASSED the 75% rule
// in red (appearing magenta over DAPI) and only FINAL accepted CAR4
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

// NEW GATE overlay: shows only nuclei that PASSED the 75% criterion,
// alongside ALL CAR4 candidates, independent of distance. If a nucleus
// appears magenta in the first candidate overlay but not in this image,
// it FAILED the nuclear-coverage gate: proximity is not the explanation.
selectWindow(ergPosMask);
run("Select None");
selectWindow(car4FilteredMask);
run("Select None");
selectWindow(dapiMask);
run("Select None");
run("Merge Channels...", "c1=[" + ergPosMask + "] c2=[" + car4FilteredMask + "] c3=[" + dapiMask + "] create keep");
rename(tripleGateQC);
run("Select None");

// NEW DISTANCE overlay: red=5-um circular neighborhood around accepted
// ERG nuclei, green=all candidate CAR4 pixels, blue=accepted ERG nuclei.
// Yellow (red+green) means a CAR4 candidate is inside the distance gate.
// Green without red means excluded by the proximity rule.
selectWindow(ergZone);
run("Select None");
selectWindow(car4FilteredMask);
run("Select None");
selectWindow(ergPosMask);
run("Select None");
run("Merge Channels...", "c1=[" + ergZone + "] c2=[" + car4FilteredMask + "] c3=[" + ergPosMask + "] create keep");
rename(tripleDistanceQC);
run("Select None");

print("NEW nuclear-gate QC (DAPI blue, accepted ERG magenta, ALL candidate CAR4 green): " + tripleGateQC);
print("NEW distance QC (red=true 5-um ERG zone, green=ALL candidate CAR4, blue=accepted ERG nuclei; yellow=within 5 um): " + tripleDistanceQC);
print("3-color QC (blue=DAPI, red=thresholded ERG inside DAPI, green=ALL filtered CAR4 candidates): " + tripleCandidateQC);
print("3-color QC (blue=DAPI, magenta=ERG-positive nuclear ROIs passing 75%, green=ACCEPTED CAR4): " + tripleAcceptedQC);
print("The new QC overlays are visualization only. IMPORTANT: this macro corrects the former square 5-pixel proximity gate to an exact circular 5-um gate, so association counts may differ from the previous macro.");

// ==================================================
// 8D. CAR4-POSITIVE ENDOTHELIAL CELL COUNT (one count per ERG+ nucleus)
// ==================================================
// This restores the CAR4 cell-count definition used in the original CAR4final(2).ijm,
// but retains the new DAPI+/ERG+ >=75% gate and the true calibrated 5-um radius.
// A CAR4-positive "cell" is operationally an ERG+ nucleus whose NON-NUCLEAR
// 5-um local region has >=5% pixels with CAR4 intensity > a distant-background
// P99 threshold. This does NOT require >=5% of the nuclear area to be CAR4+.
// Biological caveat: spatial association alone cannot assign all membrane CAR4
// signal uniquely to one endothelial cell in confluent alveolar capillary networks.
// Unlike the object-count mask, this uses the P99 threshold, not mean+3SD.
// Both output types are kept separate and clearly labeled.

nonNuclearMask = prefix + "_NONNUCLEAR_valid_mask";
car4CellBgExclusion = prefix + "_CAR4cell_background_10um_EXCLUSION";
car4CellBgValid = prefix + "_CAR4cell_distant_BACKGROUND_valid";
car4CellThresholdMask = prefix + "_CAR4cell_P99_threshold_mask";
car4CellCountNonNuc = prefix + "_CAR4cell_P99_NONNUCLEAR_mask";
car4CellPositiveMask = prefix + "_CAR4_POSITIVE_ERG_NUCLEI_5pct_P99";
car4CellQC = prefix + "_QC_DAPIblue_ERGred_CAR4P99green";
car4CellPerNucleusMask = prefix + "_TEMP_CELL_5um_zone";

// DAPI-positive regions must be excluded from the CAR4 cell-local fraction.
selectWindow(dapiMask);
run("Select None");
run("Duplicate...", "title=" + nonNuclearMask);
selectWindow(nonNuclearMask);
run("Invert");
run("Select None");

// Generate a true Euclidean 10-um exclusion region around ERG-positive nuclei
// to define distant CAR4 background, with no ROI Enlarge or Fill operations.
newImage(car4CellBgExclusion, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
bgExclusionPixels = -floor(-(car4CellBgExclusion_um / pixelSize_um)); // ImageJ macro compatible ceiling

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
        selectWindow(car4CellBgExclusion);
        for (cj=0; cj<nucX.length; cj++) {
            xx = nucX[cj]; yy = nucY[cj];
            for (ddy=-bgExclusionPixels; ddy<=bgExclusionPixels; ddy++) {
                qy = yy + ddy;
                if (qy>=0 && qy<imgHeight) {
                    for (ddx=-bgExclusionPixels; ddx<=bgExclusionPixels; ddx++) {
                        qx = xx + ddx;
                        r2 = (ddx*pixelWidth_um)*(ddx*pixelWidth_um) + (ddy*pixelHeight_um)*(ddy*pixelHeight_um);
                        if (qx>=0 && qx<imgWidth && r2<=car4CellBgExclusion_um*car4CellBgExclusion_um)
                            setPixel(qx,qy,255);
                    }
                }
            }
        }
    }
}
selectWindow(car4CellBgExclusion);
run("Select None");
run("Invert"); // White = outside 10-um exclusion zone
selectWindow(nonNuclearMask);
run("Select None");
imageCalculator("AND create", nonNuclearMask, car4CellBgExclusion);
rename(car4CellBgValid);
run("Select None");

// Measure the histogram ONLY in the valid distant-background pixels.
// Crucially, this is a full 16-bit histogram (65,536 bins), not the legacy
// 256-bin, 8-bit histogram. With 16-bit data, bin index = intensity value.
selectWindow(car4CellBgValid);
run("Select None");
getHistogram(validBgVals, validBgCounts, 256);
expectedBgPixels = validBgCounts[255];
run("Create Selection");
if (selectionType() == -1)
    exit("CELL-COUNT QC ERROR: no distant non-nuclear CAR4 background available.");
roiManager("Reset");
roiManager("Add");
selectWindow(car4Bgsub);
roiManager("Select", 0);
if (imageBits == 16) {
    getHistogram(histValues, histCounts, 65536);
    nHistBins = 65536;
} else {
    getHistogram(histValues, histCounts, 256);
    nHistBins = 256;
}
roiManager("Reset");
selectWindow(car4Bgsub);
run("Select None");
bgTotal = 0;
for (hi=0; hi<nHistBins; hi++) bgTotal += histCounts[hi];
if (bgTotal < 100)
    exit("CELL-COUNT QC ERROR: insufficient distant background pixels (<100).");
if (abs(bgTotal - expectedBgPixels) > 1)
    exit("CELL-COUNT QC ERROR: background selection included " + bgTotal + " pixels, but valid-background mask has " + expectedBgPixels + ". Check mask ROI fidelity.");
bgCum = 0;
car4CellBgP99 = 0;
for (hi=0; hi<nHistBins; hi++) {
    bgCum += histCounts[hi];
    if (bgCum >= 0.99 * bgTotal) {
        car4CellBgP99 = hi;
        break;
    }
}
print("CAR4 cell-count distant non-nuclear background: " + bgTotal + " pixels");
print("CAR4 cell-count distant-background P99 threshold (bg-subtracted intensity): " + car4CellBgP99);

// Binary CAR4-positive pixels for the cell-count rule.
// Distinct from the mean+3SD particle threshold used for CAR4 area/object counts.
selectWindow(car4Bgsub);
run("Select None");
run("Duplicate...", "title=" + car4CellThresholdMask);
selectWindow(car4CellThresholdMask);
getStatistics(car4A, car4Mean, car4Min, car4Max);
if (car4CellBgP99 >= car4Max) {
    run("Select All"); run("Clear"); run("Select None");
    run("8-bit");
    print("WARNING: no CAR4 signal is above distant-background P99.");
} else {
    setThreshold(car4CellBgP99 + 1, car4Max);
    run("Convert to Mask");
    getStatistics(maskA, maskMean);
    if (maskMean > 127) run("Invert");
}
selectWindow(car4CellThresholdMask); run("Select None");
selectWindow(nonNuclearMask); run("Select None");
imageCalculator("AND create", car4CellThresholdMask, nonNuclearMask);
rename(car4CellCountNonNuc);
run("Select None");

// Build a TRUE 5-um neighborhood for EACH nucleus. Use a temporary binary mask
// to prevent double-counting pixels reached from multiple pixels of one nucleus.
// Count local binary signal only outside the DAPI mask.
newImage(car4CellPositiveMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
newImage(car4CellPerNucleusMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");
roiManager("Reset");
if (nERGpos > 0) roiManager("Open", ergRoiZip);
if (roiManager("count") != nERGpos)
    exit("CELL-COUNT QC ERROR: expected " + nERGpos + " ERG nuclei, but restored " + roiManager("count") + " ROIs.");

nCAR4positiveCells = 0;
nCAR4cellNoValidLocalArea = 0;
for (ci=0; ci<nERGpos; ci++) {
    selectWindow(ergPosMask);
    roiManager("Select", ci);
    if (selectionType() == -1)
        exit("CELL-COUNT QC ERROR: invalid ERG nucleus index " + ci);
    Roi.getContainedPoints(nucX, nucY);

    selectWindow(car4CellPerNucleusMask);
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
    selectWindow(car4CellCountNonNuc);
    roiManager("Select", localROIindex);
    getHistogram(posVals, posCounts, 256);
    localCAR4 = posCounts[255];
    roiManager("Select", localROIindex);
    roiManager("Delete");

    localCAR4Fraction = 0;
    if (localValid > 0) localCAR4Fraction = localCAR4 / localValid;
    else nCAR4cellNoValidLocalArea++;
    if (localValid > 0 && localCAR4Fraction >= car4CellMinLocalFraction) {
        nCAR4positiveCells++;
        selectWindow(car4CellPositiveMask);
        for (cj=0; cj<nucX.length; cj++)
            setPixel(nucX[cj], nucY[cj], 255);
    }
}
selectWindow(car4CellPositiveMask); run("Select None");
// Keep car4CellPerNucleusMask open for the DAPI+/CAR4+ gate below.
if (nCAR4cellNoValidLocalArea > 0)
    print("QC WARNING: " + nCAR4cellNoValidLocalArea + " nuclei had no valid non-nuclear pixels in the local region.");
if (nCAR4positiveCells > nERGpos)
    exit("CELL-COUNT QC ERROR: CAR4-positive cells exceed ERG-positive nuclei.");
car4PositiveCellFraction = 0;
if (nERGpos > 0) car4PositiveCellFraction = nCAR4positiveCells/nERGpos;

// Image for QC: blue DAPI, magenta DAPI+/ERG+ nuclei, green CAR4 above P99
// in non-nuclear pixels. The separate CAR4-positive-nuclei mask shows which
// nuclei passed the >=5% gate.
selectWindow(ergPosMask); run("Select None");
selectWindow(car4CellCountNonNuc); run("Select None");
selectWindow(dapiMask); run("Select None");
run("Merge Channels...", "c1=[" + ergPosMask + "] c2=[" + car4CellCountNonNuc + "] c3=[" + dapiMask + "] create keep");
rename(car4CellQC);
run("Select None");

// ==================================================
// 8E. DAPI+/CAR4+ nuclei irrespective of ERG status
// ==================================================
// Apply EXACTLY the same CAR4 P99 threshold, 5-um local circular region,
// non-nuclear valid-pixel denominator, and >=5% rule as in section 8D.
// Nucleus is counted once, independently of the number of CAR4 fragments.
// The comparison is paired: every ERG+/CAR4+ nucleus must also appear here.
// This is an operational spatial-association count, not definite cell membrane
// assignment or proof that a CAR4-positive cell is ERG negative.

dapiCAR4posMask = prefix + "_DAPI_CAR4_POSITIVE_NUCLEI_5pct_P99";
dapiCAR4QC = prefix + "_QC_DAPIblue_CAR4posNUCLEIred_CAR4P99green";
newImage(dapiCAR4posMask, "8-bit black", imgWidth, imgHeight, 1);
run("Set Scale...", "distance=1 known=" + pixelSize_um + " unit=micron global");

roiManager("Reset");
roiManager("Open", dapiRoiZip);
if (roiManager("count") != nDAPI)
    exit("DAPI/CAR4 CELL-COUNT QC ERROR: expected " + nDAPI + " DAPI nuclear ROIs, restored " + roiManager("count") + ".");

nDAPICAR4positiveCells = 0;
nDAPICAR4noValidLocalArea = 0;
for (ci=0; ci<nDAPI; ci++) {
    selectWindow(dapiMask);
    roiManager("Select", ci);
    if (selectionType() == -1)
        exit("DAPI/CAR4 QC ERROR: invalid DAPI nucleus at ROI index " + ci);
    Roi.getContainedPoints(nucX, nucY);
    if (nucX.length == 0)
        exit("DAPI/CAR4 QC ERROR: empty DAPI nucleus at ROI index " + ci);

    selectWindow(car4CellPerNucleusMask);
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
        exit("DAPI/CAR4 QC ERROR: no local zone for DAPI ROI index " + ci);
    roiManager("Add");
    localROIindex = roiManager("count") - 1;
    selectWindow(nonNuclearMask);
    roiManager("Select", localROIindex);
    getHistogram(validVals, validCounts, 256);
    localValid = validCounts[255];
    selectWindow(car4CellCountNonNuc);
    roiManager("Select", localROIindex);
    getHistogram(posVals, posCounts, 256);
    localCAR4 = posCounts[255];
    roiManager("Select", localROIindex);
    roiManager("Delete");

    localCAR4Fraction = 0;
    if (localValid > 0) localCAR4Fraction = localCAR4 / localValid;
    else nDAPICAR4noValidLocalArea++;
    if (localValid > 0 && localCAR4Fraction >= car4CellMinLocalFraction) {
        nDAPICAR4positiveCells++;
        selectWindow(dapiCAR4posMask);
        for (cj=0; cj<nucX.length; cj++)
            setPixel(nucX[cj], nucY[cj], 255);
    }
}
selectWindow(dapiCAR4posMask); run("Select None");
selectWindow(car4CellPerNucleusMask); close();

if (nDAPICAR4noValidLocalArea > 0)
    print("QC WARNING: " + nDAPICAR4noValidLocalArea + " DAPI nuclei lacked a non-nuclear local region for CAR4 testing.");
if (nDAPICAR4positiveCells > nDAPI)
    exit("DAPI/CAR4 CELL-COUNT QC ERROR: positive count exceeds number of DAPI nuclei.");
if (nCAR4positiveCells > nDAPICAR4positiveCells)
    exit("DAPI/CAR4 CELL-COUNT QC ERROR: ERG+/CAR4+ cell count exceeds total DAPI+/CAR4+ count, despite identical rules. Check ROI and mask alignment.");

dapiCAR4fraction = nDAPICAR4positiveCells / nDAPI;
nDAPICAR4_ERGnegative = nDAPICAR4positiveCells - nCAR4positiveCells;

selectWindow(dapiCAR4posMask); run("Select None");
selectWindow(car4CellCountNonNuc); run("Select None");
selectWindow(dapiMask); run("Select None");
run("Merge Channels...", "c1=[" + dapiCAR4posMask + "] c2=[" + car4CellCountNonNuc + "] c3=[" + dapiMask + "] create keep");
rename(dapiCAR4QC);
run("Select None");

File.delete(dapiRoiZip);
File.delete(ergRoiZip);
print("DAPI+/CAR4+ NUCLEI, regardless of ERG (5% local/P99): " + nDAPICAR4positiveCells + " / " + nDAPI + " = " + dapiCAR4fraction);
print("DAPI+/CAR4+ nuclei that FAILED the ERG gate: " + nDAPICAR4_ERGnegative + " (not proven ERG-negative cells).");
print("New DAPI/CAR4 cell-count QC images: " + dapiCAR4posMask + " and " + dapiCAR4QC);
print("CAR4-POSITIVE ERG+ ENDOTHELIAL CELLS (cell-level 5%/P99 rule): " + nCAR4positiveCells + " / " + nERGpos + " = " + car4PositiveCellFraction);
print("Cell rule: >= " + 100*car4CellMinLocalFraction + "% of non-nuclear pixels within circular " + nearDistance_um + " um region > same-image distant-background P99.");
print("These are cell-level CAR4-associated estimates, NOT counts of CAR4 mask fragments.");
print("Cell count QC images: " + car4CellPositiveMask + " and " + car4CellQC);

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
print("CAR4-proximity-associated ERG-positive nuclei: " + nCAR4assocERG);
print("CAR4-proximity-associated fraction of ERG-positive nuclei: " + car4AssocERGFraction);
print("DAPI+/CAR4+ nuclei regardless of ERG (5% local / P99): " + nDAPICAR4positiveCells);
print("DAPI+/CAR4+ fraction of all DAPI nuclei (5% local / P99): " + dapiCAR4fraction);
print("DAPI+/CAR4+ nuclei without >=75% ERG overlap: " + nDAPICAR4_ERGnegative);
print("DAPI+/ERG+/CAR4+ nuclei (5% local / P99): " + nCAR4positiveCells);
print("CAR4-positive fraction of ERG-positive nuclei (5% local / P99): " + car4PositiveCellFraction);
print("CAR4-positive cell rule: >=5% non-nuclear pixels within true 5 um > distant-background P99=" + car4CellBgP99);
print("CAR4 threshold method: " + car4ThresholdMethod);
print("CAR4 candidate objects before ERG-distance filtering: " + nCAR4candidate);
print("ERG-distance-associated CAR4 objects: " + nAccepted);
print("Total ERG-distance-associated CAR4 area, um2: " + totalArea);
print("Total raw CAR4 integrated density: " + totalRawIntegratedDensity);
print("Total background-subtracted CAR4 integrated density: " + totalBgSubIntegratedDensity);

if (nERGpos > 0) {
    print("ERG-distance-associated CAR4 objects per ERG-positive nucleus: " + nAccepted / nERGpos);
    print("Background-subtracted CAR4 integrated density per ERG-positive nucleus: " + totalBgSubIntegratedDensity / nERGpos);
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
qcDAPIName = "QC1_CAR4_DAPI_BLUE_n" + nDAPI + "_" + qcRunTag;
qcERGName = "QC2_CAR4_DAPIblue_ERGmagenta_n" + nERGpos + "_" + qcRunTag;
qcTripleName = "QC3_CAR4_TRIPLEmagenta_CAR4green_n" + nCAR4positiveCells + "_" + qcRunTag;

// Verify segmentation masks and above-P99 marker map are available.
qcTitlesBefore = getList("image.titles");
qcDAPIok = 0; qcERGok = 0; qcTripleOK = 0; qcMarkerOK = 0;
for (qi=0; qi<qcTitlesBefore.length; qi++) {
    if (qcTitlesBefore[qi] == dapiDetectedQCMask) qcDAPIok = 1;
    if (qcTitlesBefore[qi] == ergPosMask) qcERGok = 1;
    if (qcTitlesBefore[qi] == car4CellPositiveMask) qcTripleOK = 1;
    if (qcTitlesBefore[qi] == car4CellCountNonNuc) qcMarkerOK = 1;
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
// CAR4 pixels above background P99 are green as CONTEXT, not cell counts.
// Only magenta nuclear regions in QC3 count towards the triple-positive total.
selectWindow(car4CellPositiveMask); run("Select None");
selectWindow(car4CellCountNonNuc); run("Select None");
selectWindow(dapiDetectedQCMask); run("Select None");
run("Merge Channels...", "c1=[" + car4CellPositiveMask + "] c2=[" + car4CellCountNonNuc + "] c3=[" + dapiDetectedQCMask + "] create keep");
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
print("QC3: Blue DAPI + magenta triple-positive nuclei n=" + nCAR4positiveCells + " + green above-P99 CAR4 context (" + qcTripleName + ")");
print("Closed " + closedQC + " intermediate images. Original TIFF remains open.");
print("GREEN represents CAR4 fluorescence above P99, NOT separate cells; only MAGENTA nuclei in QC3 are counted as triple-positive.");
