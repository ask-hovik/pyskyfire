FeatureScript 3070;
import(path : "onshape/std/common.fs", version : "3070.0");

// Schema v8: plain hot-wall spans revolve a spline fitted through the chamber
// contour, one per smooth run of it. Rounded and mixed spans loft one curve
// per channel/blend along independent profile runs.
// Uncovered gaps and margins around interleave boundaries remain open.
// Geometry is resolved in Python; this feature only constructs the export.

function arcPoint(arc is map, fraction is number)
{
    if (arc.kind == "line")
    {
        return ((1 - fraction) * vector(arc.start_uv_mm[0], arc.start_uv_mm[1]) +
                fraction * vector(arc.end_uv_mm[0], arc.end_uv_mm[1])) * millimeter;
    }
    const angle = (arc.startAngle_deg + fraction * (arc.endAngle_deg - arc.startAngle_deg)) * degree;
    return (vector(arc.center_uv_mm[0], arc.center_uv_mm[1]) +
            cos(angle) * vector(arc.axisA_uv_mm[0], arc.axisA_uv_mm[1]) +
            sin(angle) * vector(arc.axisB_uv_mm[0], arc.axisB_uv_mm[1])) * millimeter;
}

function addArc(sketch is Sketch, arcId is string, arc is map)
{
    if (arc.kind == "line")
    {
        skLineSegment(sketch, arcId, { "start" : arcPoint(arc, 0), "end" : arcPoint(arc, 1) });
    }
    else if (arc.kind == "ellipse")
    {
        skEllipticalArc(sketch, arcId, {
            "center" : vector(arc.center_uv_mm[0], arc.center_uv_mm[1]) * millimeter,
            "majorAxis" : vector(arc.majorAxis[0], arc.majorAxis[1]),
            "majorRadius" : arc.majorRadius_mm * millimeter,
            "minorRadius" : arc.minorRadius_mm * millimeter,
            "startParameter" : arc.startParameter,
            "endParameter" : arc.endParameter
        });
    }
    else if (arc.kind == "circle")
    {
        skArc(sketch, arcId, {
            "start" : arcPoint(arc, 0),
            "mid" : arcPoint(arc, 0.5),
            "end" : arcPoint(arc, 1)
        });
    }
    else
    {
        throw regenError("Unsupported arc kind");
    }
}

function addHotWallCurve(sketch is Sketch, curveId is string, curve is map)
{
    if (curve.kind == "spline")
    {
        // One interpolated spline per smooth contour run, so the revolved
        // wall is a single smooth face rather than a chain of segments.
        var points = [];
        for (var point in curve.points_uv_mm)
        {
            points = append(points, vector(point[0], point[1]) * millimeter);
        }
        skFitSpline(sketch, curveId, { "points" : points });
    }
    else if (curve.kind == "bezier")
    {
        var points = [];
        for (var point in curve.controlPoints_uv_mm)
        {
            points = append(points, vector(point[0], point[1]) * millimeter);
        }
        skBezier(sketch, curveId, { "points" : points });
    }
    else
    {
        addArc(sketch, curveId, curve);
    }
}

// Query the solved region, not every edge created by the sketch. Sketches
// contain both wire geometry and region topology; the loft needs one profile.
function checkedSketchRegion(context is Context, sketchId is Id, label is string)
{
    const regions = qSketchRegion(sketchId);
    const count = size(evaluateQuery(context, regions));
    if (count != 1)
    {
        throw regenError(label ~ " must enclose exactly one region; found " ~ count ~
                         ". Check the section curves for gaps or intersections.");
    }
    return regions;
}

// One sketch curve as one wire edge. A sketch that closes a region also owns
// face edges along the same geometry, so restrict this to the wire body.
function hotWallCurveEdge(context is Context, sketchId is Id, curveIndex is number, label is string)
{
    const edges = qBodyType(
        sketchEntityQuery(sketchId, EntityType.EDGE, "curve" ~ curveIndex),
        BodyType.WIRE
    );
    const count = size(evaluateQuery(context, edges));
    if (count != 1)
    {
        throw regenError(label ~ " resolved to " ~ count ~ " edges instead of one; " ~
                         "the hot-wall sketch did not create one edge per exported curve.");
    }
    return edges;
}

function hotWallSectionSketch(context is Context, sketchId is Id, section is map)
{
    const origin = section.plane.origin_mm;
    var sketch = newSketchOnPlane(context, sketchId, {
        "sketchPlane" : plane(vector(origin[0], origin[1], origin[2]) * millimeter,
                              vector(0, 0, 1), vector(1, 0, 0))
    });
    for (var j = 0; j < size(section.curves); j += 1)
    {
        addHotWallCurve(sketch, "curve" ~ j, section.curves[j]);
    }
    skSolve(sketch);
}

// Loft one exported curve along every section of one run. Each profile is a
// single edge, so there are no interior vertices left for opLoft to pair up.
function loftHotWallStrip(
    context is Context,
    stripId is Id,
    sketchIds is array,
    sectionIndices is array,
    curveIndex is number,
    curveName is string)
{
    var profiles = [];
    for (var k = 0; k < size(sectionIndices); k += 1)
    {
        const sectionIndex = sectionIndices[k];
        profiles = append(profiles, hotWallCurveEdge(
            context, sketchIds[sectionIndex], curveIndex,
            "Hot-wall curve " ~ curveName ~ " on section " ~ (sectionIndex + 1)
        ));
    }
    try
    {
        opLoft(context, stripId, {
            "profileSubqueries" : profiles,
            "bodyType" : ToolBodyType.SURFACE
        });
    }
    catch (error)
    {
        println("PySkyFire hot-wall strip loft failed for curve " ~ curveName ~ ":");
        println(error);
        throw regenError("Hot-wall loft failed for curve " ~ curveName ~
                         " over " ~ size(sectionIndices) ~ " sections. " ~
                         "See FeatureScript notices for the kernel error.");
    }
}

function knitHotWall(context is Context, id is Id, surfaces is array, knit is boolean)
{
    if (knit && size(surfaces) > 1)
    {
        try
        {
            opBoolean(context, id + "knit", {
                "tools" : qUnion(surfaces),
                "operationType" : BooleanOperationType.UNION,
                "allowSheets" : true
            });
        }
        catch (error)
        {
            println("PySkyFire could not knit this hot-wall span:");
            println(error);
            reportFeatureInfo(context, id, "Hot-wall span built as separate surfaces; knitting failed.");
        }
    }
}

function revolveHotWall(context is Context, id is Id, profile is map)
{
    const sketchId = id + "profile";
    // Sketch u = radius (world X), v = axial (world Z).
    var sketch = newSketchOnPlane(context, sketchId, {
        "sketchPlane" : plane(vector(0, 0, 0) * millimeter,
                              vector(0, -1, 0), vector(1, 0, 0))
    });
    var edges = [];
    for (var j = 0; j < size(profile.curves); j += 1)
    {
        addHotWallCurve(sketch, "curve" ~ j, profile.curves[j]);
        edges = append(edges, sketchEntityQuery(sketchId, EntityType.EDGE, "curve" ~ j));
    }
    skSolve(sketch);
    opRevolve(context, id + "revolve", {
        "entities" : qUnion(edges),
        "axis" : line(vector(0, 0, 0) * millimeter, vector(0, 0, 1)),
        "angleForward" : 360 * degree
    });
}

function createHotWall(context is Context, id is Id, wall is map,
                       knit is boolean, fullContour is boolean)
{
    if (fullContour)
    {
        if (wall.fullContourAllowed != true || wall.fullContour == undefined)
        {
            throw regenError("Full-contour hot wall requires squared channels only");
        }
        revolveHotWall(context, id + "fullContour", wall.fullContour);
        return;
    }
    for (var r = 0; r < size(wall.revolves); r += 1)
    {
        revolveHotWall(context, id + ("plain" ~ r), wall.revolves[r]);
    }
    var sketchIds = [];
    for (var i = 0; i < size(wall.sections); i += 1)
    {
        const sketchId = id + ("section" ~ i);
        hotWallSectionSketch(context, sketchId, wall.sections[i]);
        checkedSketchRegion(context, sketchId, "Hot-wall section " ~ (i + 1));
        sketchIds = append(sketchIds, sketchId);
    }
    for (var r = 0; r < size(wall.runs); r += 1)
    {
        const run = wall.runs[r];
        if (size(run.sectionIndices) < 2)
        {
            throw regenError("A hot-wall loft span needs at least two profiles");
        }
        var surfaces = [];
        const runId = id + ("run" ~ r);
        for (var j = 0; j < size(run.curveNames); j += 1)
        {
            const stripId = runId + ("strip" ~ j);
            loftHotWallStrip(context, stripId, sketchIds, run.sectionIndices, j,
                             run.curveNames[j]);
            surfaces = append(surfaces, qCreatedBy(stripId, EntityType.BODY));
        }
        // Knit each continuous span independently; never connect across gaps.
        knitHotWall(context, runId, surfaces, knit);
    }
    if (size(wall.revolves) == 0 && size(wall.runs) == 0)
    {
        reportFeatureInfo(context, id, "No hot-wall spans remain outside the interleave gaps.");
    }
}

function addArcAndSideProfile(sketch is Sketch, section is map, includeSides is boolean)
{
    addArc(sketch, "innerArc", section.innerArc);
    addArc(sketch, "outerArc", section.outerArc);
    if (includeSides)
    {
        skLineSegment(sketch, "sidePlus", {
            "start" : arcPoint(section.innerArc, 1),
            "end" : arcPoint(section.outerArc, 0)
        });
        skLineSegment(sketch, "sideMinus", {
            "start" : arcPoint(section.outerArc, 1),
            "end" : arcPoint(section.innerArc, 0)
        });
    }
}

function createSectionSketch(
    context is Context,
    sketchId is Id,
    section is map,
    sectionType is string)
{
    const originData = section.plane.origin_mm;
    const normalData = section.plane.normal;
    const xAxisData = section.plane.xAxis;
    const origin = vector(originData[0], originData[1], originData[2]) * millimeter;
    const normal = normalize(vector(normalData[0], normalData[1], normalData[2]));
    const suppliedXAxis = vector(xAxisData[0], xAxisData[1], xAxisData[2]);
    const xAxis = normalize(suppliedXAxis - dot(suppliedXAxis, normal) * normal);
    const sketchPlane = plane(origin, normal, xAxis);
    var sketch = newSketchOnPlane(context, sketchId, { "sketchPlane" : sketchPlane });

    if (sectionType == "squared")
    {
        addArcAndSideProfile(sketch, section, true);
    }
    else if (sectionType == "rounded")
    {
        // At minimum rounded height the lobes form a circle. Keep the two
        // arcs for correspondence but never create zero-length side lines.
        addArcAndSideProfile(sketch, section, !section.circular);
    }
    else
    {
        throw regenError("Unsupported cooling-channel section type");
    }
    skSolve(sketch);
    return sketchPlane;
}


function patternSeedBody(
    context is Context,
    patternId is Id,
    seedBody is Query,
    patternData is map,
    patternAxis is Line)
{
    const count = patternData.count;
    if (count <= 1)
    {
        return;
    }

    var transforms = [];
    var instanceNames = [];
    for (var instance = 1; instance < count; instance += 1)
    {
        transforms = append(
            transforms,
            rotationAround(
                patternAxis,
                instance * patternData.angularStep_deg * degree
            )
        );
        instanceNames = append(instanceNames, "instance" ~ instance);
    }

    opPattern(
        context,
        patternId,
        {
            "entities" : seedBody,
            "transforms" : transforms,
            "instanceNames" : instanceNames
        }
    );
}


// Conditional UI must depend on a boolean parameter, so editing logic reads
// JSON eligibility into hidden state. Regeneration checks the JSON again.
// https://cad.onshape.com/FsDoc/uispec.html#conditional-visibility
export function hotWallEditLogic(context is Context, id is Id, oldDefinition is map,
                                 definition is map, isCreating is boolean) returns map
{
    definition.fullContourAvailable = false;
    if (definition.channelData is JSONData)
    {
        const data = definition.channelData.jsonData;
        if (data.schemaVersion == 8 && data.hotWall != undefined)
        {
            definition.fullContourAvailable = data.hotWall.fullContourAllowed == true;
        }
    }
    if (!definition.fullContourAvailable)
    {
        definition.coverEntireHotWall = false;
    }
    return definition;
}

annotation {
    "Feature Type Name" : "PySkyFire cooling channels",
    "Editing Logic Function" : "hotWallEditLogic"
}
export const pyskyfireCoolingChannels = defineFeature(
    function(context is Context, id is Id, definition is map)
    precondition
    {
        annotation { "Name" : "Channel data" }
        definition.channelData is JSONData;

        annotation {
            "Name" : "Show section planes",
            "Default" : false
        }
        definition.showPlanes is boolean;

        annotation { "Name" : "Construct hot wall", "Default" : false }
        definition.constructHotWall is boolean;

        annotation { "Name" : "Full contour available", "Default" : false,
                     "UIHint" : UIHint.ALWAYS_HIDDEN }
        definition.fullContourAvailable is boolean;

        if (definition.constructHotWall)
        {
            if (definition.fullContourAvailable)
            {
                annotation { "Name" : "Cover entire hot wall", "Default" : false }
                definition.coverEntireHotWall is boolean;
            }
            annotation { "Name" : "Knit each hot-wall span", "Default" : true }
            definition.knitHotWall is boolean;
        }
    }
    {
        const data = definition.channelData.jsonData;
        if (data.schema != "pyskyfire.onshape.cooling-channels")
        {
            throw regenError("Unsupported cooling-channel JSON schema");
        }
        if (data.schemaVersion != 8)
        {
            throw regenError("Unsupported cooling-channel JSON version; regenerate using export_onshape.py");
        }

        if (definition.constructHotWall && data.hotWall == undefined)
        {
            throw regenError("JSON has no hot-wall data; export again with --hot-wall");
        }

        const axisData = data.coordinateSystem.patternAxis;
        const axisOriginData = axisData.origin_mm;
        const axisDirectionData = axisData.direction;
        const axisOrigin = vector(
            axisOriginData[0],
            axisOriginData[1],
            axisOriginData[2]
        ) * millimeter;
        const axisDirection = normalize(vector(
            axisDirectionData[0],
            axisDirectionData[1],
            axisDirectionData[2]
        ));
        const patternAxis = line(axisOrigin, axisDirection);

        for (var channelIndex = 0;
             channelIndex < size(data.channels);
             channelIndex += 1)
        {
            const channel = data.channels[channelIndex];
            if (size(channel.sections) < 2)
            {
                throw regenError("A cooling-channel loft needs at least two sections");
            }
            var profiles = [];

            for (var sectionIndex = 0;
                 sectionIndex < size(channel.sections);
                 sectionIndex += 1)
            {
                const section = channel.sections[sectionIndex];
                const sketchId = id + (
                    "channel" ~ channelIndex ~ "Section" ~ sectionIndex
                );
                const sectionPlane = createSectionSketch(
                    context,
                    sketchId,
                    section,
                    channel.sectionType
                );

                if (definition.showPlanes)
                {
                    opPlane(
                        context,
                        id + (
                            "channel" ~ channelIndex ~ "Plane" ~ sectionIndex
                        ),
                        {
                            "plane" : sectionPlane,
                            "width" : 120 * millimeter,
                            "height" : 120 * millimeter
                        }
                    );
                }

                profiles = append(profiles, checkedSketchRegion(
                    context, sketchId, channel.name ~ " coolant section " ~ (sectionIndex + 1)
                ));
            }

            // One channel is two or four curves in a frame that follows its
            // own helix, so proximity matching has nothing to alias onto and
            // this stays a single solid loft.
            const loftId = id + ("channelLoft" ~ channelIndex);
            try
            {
                opLoft(
                    context,
                    loftId,
                    { "profileSubqueries" : profiles }
                );
            }
            catch (error)
            {
                println("PySkyFire coolant loft failed for " ~ channel.name ~ ":");
                println(error);
                throw regenError("Coolant loft failed for " ~ channel.name ~
                                 "; all section sketches enclosed one region. See FeatureScript notices for the kernel error.");
            }

            patternSeedBody(
                context,
                id + ("channelPattern" ~ channelIndex),
                qCreatedBy(loftId, EntityType.BODY),
                channel.pattern,
                patternAxis
            );
        }
        if (definition.constructHotWall)
        {
            createHotWall(context, id + "hotWall", data.hotWall,
                          definition.knitHotWall != false,
                          data.hotWall.fullContourAllowed == true && definition.coverEntireHotWall == true);
        }
    }
);
