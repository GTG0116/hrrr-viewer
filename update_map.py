<link href="https://fonts.googleapis.com/css2?family=Montserrat:wght@800;900&display=swap" rel="stylesheet">

<div class="dashboard-container">
    <div class="glass-card">
        <div class="header">
            <div class="title-group">
                <span class="sub-title">HRRR MODEL • NORTHEAST SECTOR</span>
                <h1 id="mode-display">Temperature (°F)</h1>
            </div>
            
            <select id="modeSelect" onchange="setMode(this.value)">
                <option value="frames_temp">Temperature</option>
                <option value="frames_precip">Precip Type</option>
                <option value="frames_wind">Wind Speed</option>
                <option value="frames_gust">Wind Gusts</option>
                <option value="frames_snow">Total Snowfall</option>
                <option value="frames_total_precip">Total Precip</option>
            </select>
        </div>

        <div class="viewer-area">
            <div id="loading" class="spinner"></div>
            <img id="viewer" src="https://raw.githubusercontent.com/GTG0116/hrrr-viewer/main/frames_temp/f01.png" onload="hideLoad()">
        </div>

        <div class="controls-grid">
            <div class="play-section">
                <button id="playBtn" onclick="togglePlay()" class="btn-main">
                    <span id="playIcon">▶</span> PLAY
                </button>
                <div class="frame-counter">HOUR: <span id="frame-num">f01</span></div>
            </div>

            <div class="slider-group">
                <label class="control-label">FORECAST TIMELINE</label>
                <input type="range" id="frameSlider" min="1" max="18" value="1" oninput="updateFrame(this.value)">
                <div class="slider-labels">
                    <span>Current</span>
                    <span>+18h Forecast</span>
                </div>
            </div>

            <div class="slider-group">
                <label class="control-label">LOOP SPEED</label>
                <input type="range" id="speedSlider" min="100" max="1500" value="600" step="50" dir="rtl" oninput="updateSpeed(this.value)">
                <div class="slider-labels">
                    <span>Fast</span>
                    <span>Slow</span>
                </div>
            </div>
        </div>
    </div>
</div>

<style>
    :root {
        --bg-deep: #020617;
        --navy-glass: rgba(15, 23, 42, 0.85);
        --navy-solid: #1e293b;
        --electric-blue: #38bdf8;
        --border-glass: rgba(255, 255, 255, 0.15);
        --font-heavy: 'Montserrat', sans-serif;
    }

    /* Reset & Full Background */
    body, html {
        margin: 0;
        padding: 0;
        height: 100%;
        background-color: var(--bg-deep);
    }

    .dashboard-container {
        background: radial-gradient(circle at center, #1e293b 0%, #020617 100%);
        min-height: 100vh;
        width: 100%;
        display: flex;
        justify-content: center;
        align-items: center;
        padding: 20px;
        box-sizing: border-box;
    }

    /* Modern Glass UI */
    .glass-card {
        background: var(--navy-glass);
        backdrop-filter: blur(16px);
        -webkit-backdrop-filter: blur(16px);
        border: 1px solid var(--border-glass);
        border-radius: 32px;
        width: 100%;
        max-width: 950px;
        padding: 35px;
        box-shadow: 0 40px 100px -20px rgba(0, 0, 0, 0.8);
    }

    .header {
        display: flex;
        justify-content: space-between;
        align-items: flex-end;
        margin-bottom: 30px;
    }

    .sub-title {
        color: var(--electric-blue);
        font-family: var(--font-heavy);
        font-size: 0.75rem;
        letter-spacing: 4px;
        font-weight: 900;
        text-transform: uppercase;
        margin-bottom: 6px;
        display: block;
    }

    /* THICK FONT MATCHING EPHRATA WEATHER LOGO */
    h1 {
        color: #ffffff;
        margin: 0;
        font-family: var(--font-heavy);
        font-size: 2.2rem;
        font-weight: 900;
        letter-spacing: -1.5px;
        text-transform: uppercase;
    }

    /* DROPDOWN VISIBILITY FIX */
    select {
        background: var(--navy-solid);
        color: white;
        font-family: var(--font-heavy);
        font-weight: 800;
        border: 2px solid var(--electric-blue);
        padding: 12px 20px;
        border-radius: 14px;
        outline: none;
        cursor: pointer;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }

    select option {
        background: #0f172a; /* Force dark background for options */
        color: white;
        padding: 10px;
    }

    /* Map Display */
    .viewer-area {
        position: relative;
        background: #000;
        border-radius: 24px;
        overflow: hidden;
        border: 1px solid rgba(255, 255, 255, 0.1);
        aspect-ratio: 16 / 10;
        box-shadow: inset 0 0 40px rgba(0,0,0,1);
    }

    #viewer {
        width: 100%;
        height: 100%;
        object-fit: contain;
    }

    /* Playback Grid */
    .controls-grid {
        margin-top: 35px;
        display: grid;
        grid-template-columns: 180px 1fr 200px;
        gap: 30px;
        align-items: center;
    }

    .btn-main {
        background: var(--electric-blue);
        color: #020617;
        font-family: var(--font-heavy);
        font-weight: 900;
        border: none;
        padding: 14px;
        border-radius: 16px;
        cursor: pointer;
        width: 100%;
        transition: 0.3s;
        box-shadow: 0 4px 15px rgba(56, 189, 248, 0.3);
    }

    .btn-main:hover {
        transform: scale(1.05);
        background: #7dd3fc;
    }

    .frame-counter {
        color: #64748b;
        font-family: var(--font-heavy);
        font-weight: 900;
        font-size: 0.8rem;
        text-align: center;
        margin-top: 10px;
    }

    .control-label {
        display: block;
        color: var(--electric-blue);
        font-family: var(--font-heavy);
        font-weight: 900;
        font-size: 0.65rem;
        margin-bottom: 12px;
        text-transform: uppercase;
        letter-spacing: 1px;
    }

    /* High-End Sliders */
    input[type=range] {
        -webkit-appearance: none;
        width: 100%;
        background: transparent;
    }

    input[type=range]::-webkit-slider-runnable-track {
        height: 8px;
        background: rgba(255, 255, 255, 0.1);
        border-radius: 4px;
    }

    input[type=range]::-webkit-slider-thumb {
        -webkit-appearance: none;
        height: 24px;
        width: 24px;
        border-radius: 50%;
        background: #fff;
        border: 4px solid var(--electric-blue);
        cursor: pointer;
        margin-top: -8px;
        box-shadow: 0 0 20px rgba(56, 189, 248, 0.6);
    }

    .slider-labels {
        display: flex;
        justify-content: space-between;
        color: #475569;
        font-family: var(--font-heavy);
        font-weight: 800;
        font-size: 0.7rem;
        margin-top: 12px;
    }

    .spinner {
        display: none;
        position: absolute;
        top: 50%;
        left: 50%;
        width: 60px;
        height: 60px;
        border: 6px solid rgba(255, 255, 255, 0.05);
        border-top: 6px solid var(--electric-blue);
        border-radius: 50%;
        animation: spin 0.8s linear infinite;
        transform: translate(-50%, -50%);
    }

    @keyframes spin { 100% { transform: translate(-50%, -50%) rotate(360deg); } }

    @media (max-width: 800px) {
        .controls-grid { grid-template-columns: 1fr; gap: 25px; }
        h1 { font-size: 1.6rem; }
    }
</style>

<script>
    let currentFolder = 'frames_temp';
    let currentFrame = 1;
    let isPlaying = false;
    let playInterval;
    let loopDelay = 600; // Default speed in ms

    const modeNames = {
        'frames_temp': 'Temperature (°F)',
        'frames_precip': 'Precip Type',
        'frames_wind': 'Wind Speed (kts)',
        'frames_gust': 'Wind Gusts (kts)',
        'frames_snow': 'Total Snowfall (in)',
        'frames_total_precip': 'Total Precip'
    };

    function setMode(folder) {
        currentFolder = folder;
        document.getElementById('mode-display').innerText = modeNames[folder];
        update();
    }

    function updateFrame(val) {
        currentFrame = parseInt(val);
        document.getElementById('frame-num').innerText = `f${currentFrame.toString().padStart(2, '0')}`;
        update();
    }

    function updateSpeed(val) {
        loopDelay = parseInt(val);
        // If playing, restart interval with new speed immediately
        if (isPlaying) {
            clearInterval(playInterval);
            startLoop();
        }
    }

    function update() {
        document.getElementById('loading').style.display = 'block';
        const frameStr = currentFrame.toString().padStart(2, '0');
        const url = `https://raw.githubusercontent.com/GTG0116/hrrr-viewer/main/${currentFolder}/f${frameStr}.png`;
        document.getElementById('viewer').src = url;
    }

    function hideLoad() {
        document.getElementById('loading').style.display = 'none';
    }

    function startLoop() {
        playInterval = setInterval(() => {
            currentFrame++;
            if (currentFrame > 18) currentFrame = 1;
            document.getElementById('frameSlider').value = currentFrame;
            document.getElementById('frame-num').innerText = `f${currentFrame.toString().padStart(2, '0')}`;
            update();
        }, loopDelay);
    }

    function togglePlay() {
        const btn = document.getElementById('playBtn');
        if (isPlaying) {
            clearInterval(playInterval);
            btn.innerHTML = '<span>▶</span> PLAY';
            isPlaying = false;
        } else {
            isPlaying = true;
            btn.innerHTML = '<span>⏸</span> PAUSE';
            startLoop();
        }
    }
</script>
