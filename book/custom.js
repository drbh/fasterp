// Add playground button to the menu bar
(function() {
    'use strict';

    // Wait for the DOM to be ready
    document.addEventListener('DOMContentLoaded', function() {
        // Find the right-buttons container
        var rightButtons = document.querySelector('.right-buttons');
        if (!rightButtons) return;

        // Create playground button
        var playgroundBtn = document.createElement('a');
        // Calculate path to playground from current page
        var pathToRoot = window.location.pathname.split('/').filter(x => x).length - 1;
        var prefix = pathToRoot > 0 ? '../'.repeat(pathToRoot) : './';
        playgroundBtn.href = prefix + 'playground/';
        playgroundBtn.className = 'playground-btn';
        playgroundBtn.title = 'Try in Browser';
        playgroundBtn.setAttribute('aria-label', 'Try in Browser');
        playgroundBtn.innerHTML = `
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                <polygon points="12 2 2 7 12 12 22 7 12 2"/>
                <polyline points="2 17 12 22 22 17"/>
                <polyline points="2 12 12 17 22 12"/>
            </svg>
            <span>Playground</span>
        `;

        // Insert before the print button (first child)
        rightButtons.insertBefore(playgroundBtn, rightButtons.firstChild);
    });
})();
