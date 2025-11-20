// Add custom navigation buttons to the menu bar
(function() {
    'use strict';

    // Wait for the DOM to be ready
    document.addEventListener('DOMContentLoaded', function() {
        // Find the left-buttons container
        var leftButtons = document.querySelector('.left-buttons');
        if (!leftButtons) return;

        // Create the buttons container
        var navButtons = document.createElement('div');
        navButtons.className = 'nav-buttons';

        // GitHub button
        var githubBtn = document.createElement('a');
        githubBtn.href = 'https://github.com/drbh/fasterp';
        githubBtn.target = '_blank';
        githubBtn.className = 'nav-btn';
        githubBtn.title = 'View on GitHub';
        githubBtn.innerHTML = `
            <svg viewBox="0 0 24 24" fill="currentColor">
                <path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
            </svg>
            <span>GitHub</span>
        `;

        // Playground button
        var playgroundBtn = document.createElement('a');
        // Calculate path to playground from current page
        var pathToRoot = window.location.pathname.split('/').filter(x => x).length - 1;
        var prefix = pathToRoot > 0 ? '../'.repeat(pathToRoot) : './';
        playgroundBtn.href = prefix + 'playground/';
        playgroundBtn.className = 'nav-btn';
        playgroundBtn.title = 'Try in Browser';
        playgroundBtn.innerHTML = `
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                <polygon points="12 2 2 7 12 12 22 7 12 2"/>
                <polyline points="2 17 12 22 22 17"/>
                <polyline points="2 12 12 17 22 12"/>
            </svg>
            <span>Playground</span>
        `;

        // Add buttons to container
        navButtons.appendChild(githubBtn);
        navButtons.appendChild(playgroundBtn);

        // Insert after the existing buttons
        leftButtons.appendChild(navButtons);
    });
})();
