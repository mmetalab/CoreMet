/**
 * Cytoscape node tooltip positioning — follows the mouse cursor.
 * Content and visibility are managed by the Dash server callback
 * (via mouseoverNodeData); this script handles:
 *   1. Positioning the tooltip near the cursor
 *   2. Hiding the tooltip when the mouse leaves the canvas
 */
(function () {
    "use strict";

    var TIP_OFFSET_X = 14;
    var TIP_OFFSET_Y = 14;
    var _boundWrapper = null;

    function tryAttach() {
        var wrapper = document.getElementById("cyto-canvas-wrapper");
        var tip = document.getElementById("cyto-tooltip");
        if (!wrapper || !tip) { _boundWrapper = null; return; }
        if (_boundWrapper === wrapper) return;
        _boundWrapper = wrapper;

        wrapper.addEventListener("mousemove", function (e) {
            if (tip.style.display === "none" || !tip.textContent) return;
            var rect = wrapper.getBoundingClientRect();
            var x = e.clientX - rect.left + TIP_OFFSET_X;
            var y = e.clientY - rect.top + TIP_OFFSET_Y;
            var tw = tip.offsetWidth;
            var th = tip.offsetHeight;
            if (x + tw + 4 > wrapper.offsetWidth) x = e.clientX - rect.left - tw - TIP_OFFSET_X;
            if (y + th + 4 > wrapper.offsetHeight) y = e.clientY - rect.top - th - TIP_OFFSET_Y;
            if (x < 0) x = 4;
            if (y < 0) y = 4;
            tip.style.left = x + "px";
            tip.style.top = y + "px";
        });

        // Hide tooltip when mouse leaves the canvas area
        wrapper.addEventListener("mouseleave", function () {
            tip.style.display = "none";
        });
    }

    // Poll and re-attach when DOM changes (SPA navigation)
    setInterval(tryAttach, 500);
})();
