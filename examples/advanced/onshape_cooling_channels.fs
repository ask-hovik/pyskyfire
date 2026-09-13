FeatureScript 3070;
import(path : "onshape/std/common.fs", version : "3070.0");

// Schema v6: Python resolves transverse coolant arcs and optional complete
// hot-wall profiles. No thermal geometry is calculated in FeatureScript.
//
// A hot-wall profile is one curve per channel (squared) or an arc and a
// tangent blend per channel (rounded), so a fine engine reaches a hundred
// curves per profile, all alike and evenly spaced around the axis. Asking
// opLoft to match those by proximity fails outright past roughly ninety
// curves, and a helix that advances more than half a channel pitch between
// stations matches onto the neighbouring channel, twisting the wall the wrong
// way. Python therefore emits `runs` of sections whose curves are the same
// channels in the same order, and each curve is lofted along its run on its
// own. Correspondence is then given, not guessed, at any channel count or
// helix angle. Where a circuit hand-off ends one run and starts the next,
// Python re-cuts both bounding profiles onto one set of azimuths so those,
// too, arrive already paired.

function arcPoint(arc is map, fraction is number)
{
    const angle = (arc.startAngle_deg + fraction * (arc.endAngle_deg - arc.startAngle_deg)) * degree;
    return (vector(arc.center_uv_mm[0], arc.center_uv_mm[1]) +
            cos(angle) * vector(arc.axisA_uv_mm[0], arc.axisA_uv_mm[1]) +
            sin(angle) * vector(arc.axisB_uv_mm[0], arc.axisB_uv_mm[1])) * millimeter;
}

function addArc(sketch is Sketch, arcId is string, arc is map)
{
    if (arc.kind == "ellipse")
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
    if (curve.kind == "bezier")
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

// A hand-off replaces one channel pattern with another, so its two bounding
// profiles share no curves. Python re-cuts both onto one set of azimuths, so
// they arrive here with equal curve counts spanning equal angular intervals.
// Paired spans then sit at the same azimuth, which makes matching them in
// order the closest matching as well, so several spans can share one loft.
function loftHotWallBridge(context is Context, bridgeId is Id, bridge is map)
{
    const spansPerLoft = 24;
    var sketchIds = [];
    for (var s = 0; s < size(bridge.sections); s += 1)
    {
        const sketchId = bridgeId + ("profile" ~ s);
        hotWallSectionSketch(context, sketchId, bridge.sections[s]);
        sketchIds = append(sketchIds, sketchId);
    }
    const count = size(bridge.curveNames);
    var surfaces = [];
    for (var start = 0; start < count; start += spansPerLoft)
    {
        var stop = start + spansPerLoft;
        if (stop > count)
        {
            stop = count;
        }
        var profiles = [];
        for (var s = 0; s < size(sketchIds); s += 1)
        {
            var edges = [];
            for (var j = start; j < stop; j += 1)
            {
                edges = append(edges, hotWallCurveEdge(
                    context, sketchIds[s], j,
                    "Hot-wall bridge span " ~ (j + 1) ~ " on profile " ~ (s + 1)
                ));
            }
            profiles = append(profiles, qUnion(edges));
        }
        const loftId = bridgeId + ("spans" ~ start);
        try
        {
            opLoft(context, loftId, {
                "profileSubqueries" : profiles,
                "bodyType" : ToolBodyType.SURFACE
            });
        }
        catch (error)
        {
            println("PySkyFire hot-wall transition loft failed on spans " ~
                    (start + 1) ~ " to " ~ stop ~ " of " ~ count ~ ":");
            println(error);
            throw regenError("Hot-wall transition loft failed between sections " ~
                             (bridge.beforeSection + 1) ~ " and " ~
                             (bridge.afterSection + 1) ~ ", on spans " ~
                             (start + 1) ~ " to " ~ stop ~
                             ". See FeatureScript notices for the kernel error.");
        }
        surfaces = append(surfaces, qCreatedBy(loftId, EntityType.BODY));
    }
    return surfaces;
}

function createHotWall(context is Context, id is Id, wall is map, knit is boolean)
{
    if (size(wall.sections) < 2)
    {
        throw regenError("A hot-wall loft needs at least two profiles");
    }
    if (wall.runs == undefined || wall.bridges == undefined)
    {
        throw regenError("Hot-wall data has no profile runs; regenerate using export_onshape.py");
    }

    var sketchIds = [];
    for (var i = 0; i < size(wall.sections); i += 1)
    {
        const sketchId = id + ("section" ~ i);
        hotWallSectionSketch(context, sketchId, wall.sections[i]);
        // A profile that fails to close is a Python export bug worth catching
        // before any surface work starts.
        checkedSketchRegion(context, sketchId, "Hot-wall section " ~ (i + 1));
        sketchIds = append(sketchIds, sketchId);
    }

    var surfaces = [];
    for (var r = 0; r < size(wall.runs); r += 1)
    {
        const run = wall.runs[r];
        // A lone section carries no strip of its own; the bridges on either
        // side of it cover that part of the wall.
        if (size(run.sectionIndices) < 2)
        {
            continue;
        }
        for (var j = 0; j < size(run.curveNames); j += 1)
        {
            const stripId = id + ("run" ~ r ~ "strip" ~ j);
            loftHotWallStrip(context, stripId, sketchIds, run.sectionIndices, j,
                             run.curveNames[j]);
            surfaces = append(surfaces, qCreatedBy(stripId, EntityType.BODY));
        }
    }
    for (var b = 0; b < size(wall.bridges); b += 1)
    {
        const bridged = loftHotWallBridge(context, id + ("bridge" ~ b), wall.bridges[b]);
        for (var k = 0; k < size(bridged); k += 1)
        {
            surfaces = append(surfaces, bridged[k]);
        }
    }
    if (size(surfaces) == 0)
    {
        throw regenError("Hot-wall data produced no loftable spans");
    }

    // Neighbouring strips share exported endpoints exactly, so the union
    // knits them. It is the one step here that depends on the kernel matching
    // surface edges, so keep the strips rather than fail the whole feature.
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
            println("PySkyFire could not knit the hot-wall strips into one surface:");
            println(error);
            reportFeatureInfo(context, id, "Hot wall built as " ~ size(surfaces) ~
                              " separate surfaces; knitting them failed.");
        }
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


annotation {
    "Feature Type Name" : "PySkyFire cooling channels"
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

        if (definition.constructHotWall)
        {
            annotation { "Name" : "Knit hot wall into one surface", "Default" : true }
            definition.knitHotWall is boolean;
        }
    }
    {
        const data = definition.channelData.jsonData;
        if (data.schema != "pyskyfire.onshape.cooling-channels")
        {
            throw regenError("Unsupported cooling-channel JSON schema");
        }
        if (data.schemaVersion != 6)
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
            createHotWall(context, id + "hotWall", data.hotWall, definition.knitHotWall != false);
        }
    }
);
