// グローバルスコープにチャートオブジェクトを保持
let freqChart = null;
let zPlaneChart = null;
let impulseChart = null;
let inputWaveformChart = null;
let outputWaveformChart = null;
let inputFftChart = null;
let outputFftChart = null;

let audioContext = null;
let gainNode = null;
let currentSourceNode = null;
window.isPlotting3D = false;

const initialCoefficients = {
    b0: 1.0, b1: 0.0, b2: 0.0,
    a0: 1.0, a1: 0.0, a2: 0.0
};
const initialSamplingFrequency = 48000;
const initialSignalAmplitude = 1.0;

document.addEventListener('DOMContentLoaded', () => {
    // 振幅入力も 'input' イベントでリッスンするIDリストに追加
    const coeffIds = ['b0', 'b1', 'b2', 'a1', 'a2', 'samplingFrequency', 'signalAmplitude'];
    coeffIds.forEach(id => {
        const element = document.getElementById(id);
        if (element) {
            let debounceTimer;
            element.addEventListener('input', () => {
                clearTimeout(debounceTimer);
                debounceTimer = setTimeout(() => {
                    updateAllPlots();
                }, 250);
            });
        }
    });

    const resetButton = document.getElementById('resetCoeffs');
    if (resetButton) {
        resetButton.addEventListener('click', () => {
            setInitialCoefficientsAndFs();
            updateAllPlots();
        });
    }

    const playButton = document.getElementById('playButton');
    const stopButton = document.getElementById('stopButton');
    const volumeSlider = document.getElementById('volumeSlider');
    const signalTypeSelect = document.getElementById('signalType');

    if (playButton) playButton.addEventListener('click', playCurrentSignal);
    if (stopButton) stopButton.addEventListener('click', stopAudio);
    if (volumeSlider) volumeSlider.addEventListener('input', updateVolume);
    if (signalTypeSelect) signalTypeSelect.addEventListener('change', () => {
        stopAudio();
        updateAllPlots();
    });

    setInitialCoefficientsAndFs();
    initAudio();
    updateAllPlots();
});

function initAudio() {
    try {
        audioContext = new (window.AudioContext || window.webkitAudioContext)();
        gainNode = audioContext.createGain();
        gainNode.gain.value = parseFloat(document.getElementById('volumeSlider').value);
        gainNode.connect(audioContext.destination);
    } catch (e) {
        console.error("Web Audio API is not supported in this browser", e);
        const errorDiv = document.createElement('div');
        errorDiv.textContent = "お使いのブラウザはWeb Audio APIをサポートしていません。音声機能は利用できません。";
        errorDiv.style.color = 'red';
        errorDiv.style.textAlign = 'center';
        errorDiv.style.padding = '10px';
        const audioTestSection = document.querySelector('.audio-test-section');
        if (audioTestSection) {
            const firstChild = audioTestSection.firstChild;
            if (firstChild) {
                audioTestSection.insertBefore(errorDiv, firstChild);
            } else {
                audioTestSection.appendChild(errorDiv);
            }
        }
        const audioControls = document.querySelector('.audio-controls');
        if(audioControls) audioControls.style.display = 'none';
        const audioRelatedCharts = document.querySelector('.audio-related-charts');
        if(audioRelatedCharts) audioRelatedCharts.style.display = 'none';
    }
}

function setInitialCoefficientsAndFs() {
    document.getElementById('b0').value = initialCoefficients.b0;
    document.getElementById('b1').value = initialCoefficients.b1;
    document.getElementById('b2').value = initialCoefficients.b2;
    document.getElementById('a1').value = initialCoefficients.a1;
    document.getElementById('a2').value = initialCoefficients.a2;
    document.getElementById('samplingFrequency').value = initialSamplingFrequency;
    document.getElementById('signalAmplitude').value = initialSignalAmplitude;
}

function getCoefficients() {
    const b0 = parseFloat(document.getElementById('b0').value) || 0;
    const b1 = parseFloat(document.getElementById('b1').value) || 0;
    const b2 = parseFloat(document.getElementById('b2').value) || 0;
    const a0 = 1.0;
    const a1 = parseFloat(document.getElementById('a1').value) || 0;
    const a2 = parseFloat(document.getElementById('a2').value) || 0;
    return { b0, b1, b2, a0, a1, a2 };
}

function getSamplingFrequency() {
    const fs = parseFloat(document.getElementById('samplingFrequency').value) || 48000;
    return Math.max(1, fs);
}

function getSignalAmplitude() {
    const amplitudeInput = document.getElementById('signalAmplitude');
    if (!amplitudeInput) return 1.0; // 要素が存在しない場合はデフォルト値
    const amplitude = parseFloat(amplitudeInput.value);
    if (isNaN(amplitude) || amplitude < 0 || amplitude > 1.0) {
        amplitudeInput.value = initialSignalAmplitude; // 不正な値なら初期値に戻す
        return initialSignalAmplitude;
    }
    return amplitude;
}


function calculateFrequencyResponse(coeffs, fs, numPoints = 512) {
    const { b0, b1, b2, a0, a1, a2 } = coeffs;
    const magnitudes = new Array(numPoints);
    const phases = new Array(numPoints);
    const actualFrequencies = new Array(numPoints);
    const minFreq = Math.max(1, fs * 0.0001);
    const maxFreq = fs / 2;
    const logMinFreq = Math.log10(minFreq);
    const logMaxFreq = Math.log10(maxFreq);

    for (let i = 0; i < numPoints; i++) {
        const logFreq = logMinFreq + (logMaxFreq - logMinFreq) * i / (numPoints - 1);
        actualFrequencies[i] = Math.pow(10, logFreq);
        const normalizedOmega = (actualFrequencies[i] / (fs / 2)) * Math.PI;
        const cosOmega = Math.cos(normalizedOmega);
        const sinOmega = Math.sin(normalizedOmega);
        const cos2Omega = Math.cos(2 * normalizedOmega);
        const sin2Omega = Math.sin(2 * normalizedOmega);
        const numReal = b0 + b1 * cosOmega + b2 * cos2Omega;
        const numImag = -b1 * sinOmega - b2 * sin2Omega;
        const denReal = a0 + a1 * cosOmega + a2 * cos2Omega;
        const denImag = -a1 * sinOmega - a2 * sin2Omega;
        const denMagnitudeSq = denReal * denReal + denImag * denImag;

        if (denMagnitudeSq < 1e-12) {
            magnitudes[i] = 1e6;
            phases[i] = (numImag * denReal - numReal * denImag === 0) ? 0 : Math.atan2(numImag * denReal - numReal * denImag, numReal * denReal + numImag * denImag);
        } else {
            const hReal = (numReal * denReal + numImag * denImag) / denMagnitudeSq;
            const hImag = (numImag * denReal - numReal * denImag) / denMagnitudeSq;
            magnitudes[i] = Math.sqrt(hReal * hReal + hImag * hImag);
            phases[i] = Math.atan2(hImag, hReal);
        }
    }
    return { actualFrequencies, magnitudes, phases, nyquist: fs / 2, minFreqLog: logMinFreq, maxFreqLog: logMaxFreq };
}

function solveQuadratic(a, b, c) {
    if (Math.abs(a) < 1e-9) {
        if (Math.abs(b) < 1e-9) return [];
        return [{ real: -c / b, imag: 0 }];
    }
    const discriminant = b * b - 4 * a * c;
    const twoA = 2 * a;
    const roots = [];
    if (discriminant >= 0) {
        roots.push({ real: (-b + Math.sqrt(discriminant)) / twoA, imag: 0 });
        roots.push({ real: (-b - Math.sqrt(discriminant)) / twoA, imag: 0 });
    } else {
        roots.push({ real: -b / twoA, imag: Math.sqrt(-discriminant) / twoA });
        roots.push({ real: -b / twoA, imag: -Math.sqrt(-discriminant) / twoA });
    }
    return roots;
}

function calculatePolesAndZeros(coeffs) {
    const { b0, b1, b2, a0, a1, a2 } = coeffs;
    const zeros = solveQuadratic(b0, b1, b2);
    const poles = solveQuadratic(a0, a1, a2);
    return { zeros, poles };
}

function calculateImpulseResponse(coeffs, length = 64) {
    const { b0, b1, b2, a0, a1, a2 } = coeffs;
    const y = new Array(length).fill(0);
    const x = new Array(length).fill(0);
    x[0] = 1;
    for (let n = 0; n < length; n++) {
        let val = b0 * x[n];
        if (n >= 1) val += b1 * x[n - 1];
        if (n >= 2) val += b2 * x[n - 2];
        if (n >= 1) val -= a1 * y[n - 1];
        if (n >= 2) val -= a2 * y[n - 2];
        y[n] = val / a0;
    }
    return y;
}

function calculateZPlane3DData(coeffs, gridSize = 25) {
    const { b0, b1, b2, a0, a1, a2 } = coeffs;
    const xGridValues = []; const yGridValues = []; const zMagnitudes = [];
    const range = 1.5;
    for (let i = 0; i < gridSize; i++) {
        const val = -range + (2 * range * i) / (gridSize - 1);
        xGridValues.push(val); yGridValues.push(val);
    }
    for (let j = 0; j < gridSize; j++) {
        const im = yGridValues[j]; const rowMagnitudes = [];
        for (let i = 0; i < gridSize; i++) {
            const re = xGridValues[i]; const modSqZ = re * re + im * im;
            let numReal, numImag, denReal, denImag;
            if (modSqZ < 1e-9) {
                if (Math.abs(a2) > 1e-9) { numReal = b2; numImag = 0; denReal = a2; denImag = 0; }
                else if (Math.abs(a1) > 1e-9) { numReal = b1; numImag = 0; denReal = a1; denImag = 0; }
                else if (Math.abs(a0) > 1e-9) { numReal = b0; numImag = 0; denReal = a0; denImag = 0; }
                else { rowMagnitudes.push(60); continue; }
            } else {
                const zInvRe = re / modSqZ; const zInvIm = -im / modSqZ;
                const zInvSqRe = (re * re - im * im) / (modSqZ * modSqZ);
                const zInvSqIm = (-2 * re * im) / (modSqZ * modSqZ);
                numReal = b0 + b1 * zInvRe + b2 * zInvSqRe; numImag = b1 * zInvIm + b2 * zInvSqIm;
                denReal = a0 + a1 * zInvRe + a2 * zInvSqRe; denImag = a1 * zInvIm + a2 * zInvSqIm;
            }
            const denMagSq = denReal * denReal + denImag * denImag;
            if (denMagSq < 1e-12) { rowMagnitudes.push(60); }
            else {
                const hReal = (numReal * denReal + numImag * denImag) / denMagSq;
                const hImag = (numImag * denReal - numReal * denImag) / denMagSq;
                const magnitude = Math.sqrt(hReal * hReal + hImag * hImag);
                rowMagnitudes.push(Math.min(60, Math.max(-60, 20 * Math.log10(Math.max(magnitude, 1e-6)))));
            }
        }
        zMagnitudes.push(rowMagnitudes);
    }
    const unitCircle3D_x = []; const unitCircle3D_y = []; const unitCircle3D_z = [];
    const zPlaneLevel = 0;
    for (let i = 0; i <= 360; i += 5) {
        const angle = i * Math.PI / 180;
        unitCircle3D_x.push(Math.cos(angle)); unitCircle3D_y.push(Math.sin(angle)); unitCircle3D_z.push(zPlaneLevel);
    }
    return { x: xGridValues, y: yGridValues, z: zMagnitudes, unitCircle3D: {x: unitCircle3D_x, y: unitCircle3D_y, z: unitCircle3D_z} };
}

function plotFrequencyResponse(data, chartInstance, canvasId) {
    const { actualFrequencies, magnitudes, phases, minFreqLog, maxFreqLog } = data;
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
        return chartInstance; 
    }
    const ctx = canvas.getContext('2d', { willReadFrequently: true }); 
    const magDb = magnitudes.map(m => Math.min(60, Math.max(-100, 20 * Math.log10(Math.max(m, 1e-5)))));
    const phaseDeg = phases.map(p => p * 180 / Math.PI);
    const chartDataPoints = actualFrequencies.map((freq, index) => ({ x: freq, yMag: magDb[index], yPhase: phaseDeg[index] }));

    if (chartInstance) {
        chartInstance.data.datasets[0].data = chartDataPoints.map(p => ({x: p.x, y: p.yMag}));
        chartInstance.data.datasets[1].data = chartDataPoints.map(p => ({x: p.x, y: p.yPhase}));
        chartInstance.options.scales.x.min = Math.pow(10, minFreqLog);
        chartInstance.options.scales.x.max = Math.pow(10, maxFreqLog);
        chartInstance.update('none');
        return chartInstance;
    } else {
        return new Chart(ctx, {
            type: 'line',
            data: {
                datasets: [
                    { label: '振幅特性 (dB)', data: chartDataPoints.map(p => ({x: p.x, y: p.yMag})), borderColor: 'rgb(75, 192, 192)', yAxisID: 'y-axis-mag', tension: 0.1, pointRadius: 0 },
                    { label: '位相特性 (度)', data: chartDataPoints.map(p => ({x: p.x, y: p.yPhase})), borderColor: 'rgb(255, 99, 132)', yAxisID: 'y-axis-phase', tension: 0.1, pointRadius: 0 }
                ]
            },
            options: {
                responsive: true, maintainAspectRatio: false,
                scales: {
                    x: {
                        type: 'logarithmic', title: { display: true, text: '周波数 (Hz)' },
                        min: Math.pow(10, minFreqLog), max: Math.pow(10, maxFreqLog),
                        ticks: {
                            callback: function(value) {
                                const logVal = Math.log10(value);
                                if (value === 0) return '0';
                                if (logVal === Math.floor(logVal)) {
                                     if (value >= 10000 && value % 10000 === 0) return (value / 1000) + 'k';
                                     if (value >= 1000 && value % 1000 === 0) return (value / 1000) + 'k';
                                     if (value >=1 && value <=10) return value.toFixed(0);
                                     if (value < 1) return value.toFixed(1);
                                     return value.toFixed(0);
                                }
                                const firstDigit = parseFloat(value.toExponential(0)[0]);
                                if (value > 10 && (firstDigit === 1 || firstDigit === 2 || firstDigit === 5)) {
                                     if (value >= 1000) return (value / 1000).toFixed(0) + 'k';
                                     return value.toFixed(0);
                                }
                                return null;
                            },
                            autoSkip: true, maxTicksLimit: 12,
                            major: { enabled: true }
                        }
                    },
                    'y-axis-mag': { type: 'linear', position: 'left', title: { display: true, text: '振幅 (dB)' }, min: -20, max: 20 },
                    'y-axis-phase': { type: 'linear', position: 'right', title: { display: true, text: '位相 (度)' }, min: -180, max: 180, ticks: { stepSize: 90 }, grid: { drawOnChartArea: false } }
                },
                interaction: { mode: 'index', intersect: false },
            }
        });
    }
}

function plotZPlanePolesZeros(data) {
    const { poles, zeros } = data;
    const canvas = document.getElementById('zPlaneChart');
    const stabilityStatusEl = document.getElementById('stabilityStatus');
    const poleValuesTextEl = document.getElementById('poleValuesText');
    const zeroValuesTextEl = document.getElementById('zeroValuesText');
    if(!canvas) return;
    const ctx = canvas.getContext('2d', { willReadFrequently: true });

    const poleData = poles.map(p => ({ x: p.real, y: p.imag, value: `(${p.real.toFixed(3)}, ${p.imag.toFixed(3)}j)` }));
    const zeroData = zeros.map(z => ({ x: z.real, y: z.imag, value: `(${z.real.toFixed(3)}, ${z.imag.toFixed(3)}j)` }));
    const unitCircle = [];
    for (let i = 0; i <= 360; i += 5) {
        const angle = i * Math.PI / 180;
        unitCircle.push({ x: Math.cos(angle), y: Math.sin(angle) });
    }

    let isStable = true;
    let poleTexts = [];
    poles.forEach(p => {
        const magnitude = Math.sqrt(p.real * p.real + p.imag * p.imag);
        if (magnitude >= 0.9999) {
            isStable = false;
        }
        poleTexts.push(`(${p.real.toFixed(3)}, ${p.imag.toFixed(3)}j) |mag|= ${magnitude.toFixed(3)}`);
    });

    if (stabilityStatusEl) {
        stabilityStatusEl.textContent = `安定性: ${isStable ? '安定' : '不安定'}`;
        stabilityStatusEl.className = `stability-status ${isStable ? 'stable' : 'unstable'}`;
    }
    if (poleValuesTextEl) poleValuesTextEl.textContent = poles.length > 0 ? poleTexts.join('; ') : '-';
    if (zeroValuesTextEl) zeroValuesTextEl.textContent = zeros.length > 0 ? zeroData.map(z => z.value).join('; ') : '-';

    if (zPlaneChart) {
        zPlaneChart.data.datasets[0].data = unitCircle;
        zPlaneChart.data.datasets[1].data = poleData;
        zPlaneChart.data.datasets[2].data = zeroData;
        zPlaneChart.update('none');
    } else {
        zPlaneChart = new Chart(ctx, {
            type: 'scatter',
            data: {
                datasets: [
                    { label: '単位円', data: unitCircle, borderColor: 'rgba(150, 150, 150, 0.9)', borderWidth: 2.5, pointRadius: 0, showLine: true, tension: 0.1, fill: false, order: 3 },
                    { label: '極 (Poles)', data: poleData, backgroundColor: 'rgba(255, 20, 20, 0.95)', borderColor: 'rgba(180, 0, 0, 1)', borderWidth: 1.5, pointStyle: 'crossRot', radius: 10, hoverRadius: 12, order: 1 },
                    { label: '零点 (Zeros)', data: zeroData, backgroundColor: 'rgba(0, 0, 255, 0.7)', pointStyle: 'circle', radius: 6, hoverRadius: 8, borderColor: 'rgba(0,0,255,1)', borderWidth: 1, order: 2 }
                ]
            },
            options: {
                responsive: true, maintainAspectRatio: true, // aspectRatio を有効にするために true に変更
                aspectRatio: 1, // 正円表示
                scales: {
                    x: { title: { display: true, text: '実数部 (Re)' }, min: -1.5, max: 1.5, grid: { color: '#ddd' } },
                    y: { title: { display: true, text: '虚数部 (Im)' }, min: -1.5, max: 1.5, grid: { color: '#ddd' } }
                },
                plugins: {
                    legend: { display: true },
                    tooltip: { callbacks: { label: function(context) { let label = context.dataset.label || ''; if (label) { label += ': '; } if (context.raw && context.raw.value) { label += context.raw.value; } else if (context.parsed) { label += `(${context.parsed.x.toFixed(3)}, ${context.parsed.y.toFixed(3)}j)`; } return label; } } }
                }
            }
        });
    }
}

function plotImpulseResponse(impulseData) {
    const canvas = document.getElementById('impulseResponseChart');
    if(!canvas) return;
    const ctx = canvas.getContext('2d', { willReadFrequently: true });
    const labels = impulseData.map((_, i) => i);
    if (impulseChart) {
        impulseChart.data.labels = labels;
        impulseChart.data.datasets[0].data = impulseData;
        impulseChart.update('none');
    } else {
        impulseChart = new Chart(ctx, {
            type: 'bar',
            data: {
                labels: labels,
                datasets: [{ label: 'インパルス応答 y[n]', data: impulseData, backgroundColor: 'rgba(54, 162, 235, 0.6)', borderColor: 'rgb(54, 162, 235)', borderWidth: 1 }]
            },
            options: { responsive: true, maintainAspectRatio: false, scales: { x: { title: { display: true, text: 'サンプル (n)' } }, y: { title: { display: true, text: '振幅' } } } }
        });
    }
}

async function plotZPlane3D(data) {
    const { x, y, z, unitCircle3D } = data;
    const plotlyLoader = document.getElementById('plotlyLoader');
    if (plotlyLoader) plotlyLoader.style.display = 'block';
    await new Promise(resolve => setTimeout(resolve, 50));
    const plotData = [
        { type: 'surface', x: x, y: y, z: z, colorscale: 'Viridis', contours: { z: { show: true, usecolormap: true, project: { z: true } } }, showscale: true, colorbar: {title: '振幅 (dB)'} },
        { type: 'scatter3d', mode: 'lines', x: unitCircle3D.x, y: unitCircle3D.y, z: unitCircle3D.z, line: { color: 'red', width: 4 }, name: '単位円 (z=0dB)', showlegend: true }
    ];
    const layout = {
        title: 'Z平面 3D振幅特性',
        scene: { xaxis: { title: '実数部 (Re)', range: [-1.5, 1.5] }, yaxis: { title: '虚数部 (Im)', range: [-1.5, 1.5] }, zaxis: { title: '振幅 (dB)', range: [-60, 60] }, camera: { eye: { x: 1.5, y: 1.5, z: 1.2 } }, aspectmode: 'cube' },
        autosize: true, margin: { l: 40, r: 20, b: 40, t: 60 }, legend: { x: 0.8, y: 0.9 }
    };
    try { await Plotly.newPlot('zPlane3DChart', plotData, layout, {responsive: true}); }
    catch (error) { console.error("Plotly 3D plot error:", error); const chartDiv = document.getElementById('zPlane3DChart'); if(chartDiv) chartDiv.innerHTML = "<p style='color:red; text-align:center;'>3Dグラフの描画に失敗しました。</p>"; }
    finally { if (plotlyLoader) plotlyLoader.style.display = 'none'; }
}

function generateSignalBuffer(signalTypeString, durationSec, sampleRate, amplitude = 1.0) {
    if (!audioContext) {
        console.warn("AudioContext not initialized in generateSignalBuffer. Returning null.");
        return null;
    }
    const numFrames = Math.floor(durationSec * sampleRate);
    const buffer = audioContext.createBuffer(1, numFrames, sampleRate);
    const channelData = buffer.getChannelData(0);
    let signalFreq = 0;
    if (signalTypeString.startsWith('sine_')) {
        signalFreq = parseInt(signalTypeString.split('_')[1]);
        for (let i = 0; i < numFrames; i++) {
            channelData[i] = amplitude * Math.sin(2 * Math.PI * signalFreq * i / sampleRate);
        }
    } else if (signalTypeString === 'whitenoise') {
        for (let i = 0; i < numFrames; i++) {
            channelData[i] = amplitude * (Math.random() * 2 - 1);
        }
    }
    return buffer;
}

async function playCurrentSignal() {
    if (!audioContext) {
        console.warn("AudioContext not initialized. Cannot play audio.");
        return;
    }
    if (audioContext.state === 'suspended') { await audioContext.resume(); }
    if (currentSourceNode) { stopAudio(); }

    const signalType = document.getElementById('signalType').value;
    const signalAmplitude = getSignalAmplitude();
    const fs = audioContext.sampleRate;
    const duration = 2.0;
    const inputBuffer = generateSignalBuffer(signalType, duration, fs, signalAmplitude);
    if (!inputBuffer) return;

    const coeffs = getCoefficients();

    const offlineCtx = new OfflineAudioContext(1, inputBuffer.length, fs);
    const offlineSource = offlineCtx.createBufferSource();
    offlineSource.buffer = inputBuffer;
    const offlineIIRFilter = offlineCtx.createIIRFilter(
        [coeffs.b0, coeffs.b1, coeffs.b2],
        [coeffs.a0, coeffs.a1, coeffs.a2]
    );
    offlineSource.connect(offlineIIRFilter);
    offlineIIRFilter.connect(offlineCtx.destination);
    offlineSource.start();

    try {
        const renderedBuffer = await offlineCtx.startRendering();
        const outputChannelData = renderedBuffer.getChannelData(0);
        for (let i = 0; i < outputChannelData.length; i++) {
            outputChannelData[i] = Math.max(-1.0, Math.min(1.0, outputChannelData[i]));
        }
        const clippedOutputBuffer = audioContext.createBuffer(1, renderedBuffer.length, renderedBuffer.sampleRate);
        clippedOutputBuffer.copyToChannel(outputChannelData, 0);

        await plotSignalRelatedCharts(inputBuffer, clippedOutputBuffer);

        currentSourceNode = audioContext.createBufferSource();
        currentSourceNode.buffer = clippedOutputBuffer;
        currentSourceNode.connect(gainNode);
        currentSourceNode.start();
        currentSourceNode.onended = () => { currentSourceNode = null; };
    } catch (e) {
        console.error("Error processing or playing audio: ", e);
        await plotSignalRelatedCharts(inputBuffer, null);
    }
}

function stopAudio() {
    if (currentSourceNode) {
        try { currentSourceNode.stop(); } catch(e) { /* ignore */ }
        currentSourceNode.disconnect();
        currentSourceNode = null;
    }
}

function updateVolume() {
    if (gainNode) {
        gainNode.gain.value = parseFloat(document.getElementById('volumeSlider').value);
    }
}

function calculateDft(timeDomainData, sampleRate) {
    const N = timeDomainData.length;
    if (N === 0) return { frequencies: [], magnitudes: [], phases: [] };
    const magnitudes = new Array(Math.floor(N / 2)).fill(0);
    const phases = new Array(Math.floor(N / 2)).fill(0);
    const frequencies = new Array(Math.floor(N / 2)).fill(0);

    for (let k = 0; k < N / 2; k++) {
        let realSum = 0;
        let imagSum = 0;
        for (let n = 0; n < N; n++) {
            const angle = (2 * Math.PI * k * n) / N;
            realSum += timeDomainData[n] * Math.cos(angle);
            imagSum -= timeDomainData[n] * Math.sin(angle);
        }
        magnitudes[k] = (2 / N) * Math.sqrt(realSum * realSum + imagSum * imagSum);
        if (k===0 || (k === Math.floor(N/2) -1 && N%2 === 0) ) {
             magnitudes[k] = magnitudes[k]/2;
        }
        phases[k] = Math.atan2(imagSum, realSum);
        frequencies[k] = k * sampleRate / N;
    }
    return { frequencies, magnitudes, phases };
}

function plotFftChart(chartInstance, canvasId, buffer, sampleRate, chartLabel) {
    const canvas = document.getElementById(canvasId);
     if (!canvas) {
        return chartInstance;
    }
    const ctx = canvas.getContext('2d', { willReadFrequently: true });
    const fftDefaultOptions = {
        responsive: true, maintainAspectRatio: false, animation: false,
        scales: {
            x: {
                type: 'logarithmic', title: { display: true, text: '周波数 (Hz)' },
                min: Math.max(1, sampleRate / 2048),
                max: sampleRate / 2,
                ticks: {
                    callback: function(value) {
                        const logVal = Math.log10(value); if (value === 0) return '0';
                        if (logVal === Math.floor(logVal)) { if (value >= 1000) return (value / 1000) + 'k'; return value.toFixed(0); }
                        const firstDigit = parseFloat(value.toExponential(0)[0]);
                        if (value > 10 && (firstDigit === 1 || firstDigit === 2 || firstDigit === 5)) { if (value >= 1000) return (value / 1000).toFixed(0) + 'k'; return value.toFixed(0); }
                        return null;
                    },
                    autoSkip: true, maxTicksLimit: 8, major: { enabled: true }
                }
            },
            'y-axis-mag-fft': { type: 'linear', position: 'left', title: { display: true, text: '振幅 (dB)' }, min: -100, max: 0 },
            'y-axis-phase-fft': { type: 'linear', position: 'right', title: { display: true, text: '位相 (度)' }, min: -180, max: 180, ticks: { stepSize: 90 }, grid: { drawOnChartArea: false } }
        }
    };

    if (!buffer || !(buffer instanceof AudioBuffer) || buffer.length === 0 || typeof buffer.getChannelData !== 'function') {
        if (chartInstance) {
            chartInstance.data.datasets[0].data = [];
            chartInstance.data.datasets[1].data = [];
            chartInstance.update('none');
        } else {
            const datasets = [
                { label: `振幅 (dB) - ${chartLabel}`, data: [], borderColor: 'rgb(0, 200, 0)', yAxisID: 'y-axis-mag-fft', tension: 0.1, pointRadius: 0 },
                { label: `位相 (度) - ${chartLabel}`, data: [], borderColor: 'rgb(150, 100, 255)', yAxisID: 'y-axis-phase-fft', tension: 0.1, pointRadius: 0 }
            ];
            const newChart = new Chart(ctx, { type: 'line', data: { datasets }, options: fftDefaultOptions });
            if (canvasId === 'inputFftChart') inputFftChart = newChart;
            else if (canvasId === 'outputFftChart') outputFftChart = newChart;
            return newChart;
        }
        return chartInstance;
    }

    const timeDomainData = buffer.getChannelData(0);
    const fftSize = 2048;
    const dataForFft = timeDomainData.slice(0, Math.min(fftSize, timeDomainData.length));

    if (dataForFft.length === 0) {
         if (chartInstance) { chartInstance.data.datasets[0].data = []; chartInstance.data.datasets[1].data = []; chartInstance.update('none'); }
         return chartInstance;
    }

    const { frequencies, magnitudes, phases } = calculateDft(dataForFft, sampleRate);
    const magDb = magnitudes.map(m => 20 * Math.log10(Math.max(m, 1e-7)));
    const phaseDeg = phases.map(p => p * 180 / Math.PI);
    const fftDataPoints = frequencies.map((freq, index) => ({ x: freq, yMag: magDb[index], yPhase: phaseDeg[index] }));

    if (chartInstance) {
        chartInstance.data.datasets[0].data = fftDataPoints.map(p => ({x: p.x, y: p.yMag}));
        chartInstance.data.datasets[1].data = fftDataPoints.map(p => ({x: p.x, y: p.yPhase}));
        chartInstance.options.scales.x.min = Math.max(1, sampleRate / dataForFft.length);
        chartInstance.options.scales.x.max = sampleRate / 2;
        chartInstance.update('none');
        return chartInstance;
    } else {
        const newChartOptions = JSON.parse(JSON.stringify(fftDefaultOptions));
        newChartOptions.scales.x.min = Math.max(1, sampleRate / dataForFft.length);
        const datasets = [
            { label: `振幅 (dB) - ${chartLabel}`, data: fftDataPoints.map(p => ({x: p.x, y: p.yMag})), borderColor: 'rgb(0, 200, 0)', yAxisID: 'y-axis-mag-fft', tension: 0.1, pointRadius: 0 },
            { label: `位相 (度) - ${chartLabel}`, data: fftDataPoints.map(p => ({x: p.x, y: p.yPhase})), borderColor: 'rgb(150, 100, 255)', yAxisID: 'y-axis-phase-fft', tension: 0.1, pointRadius: 0 }
        ];
        return new Chart(ctx, { type: 'line', data: { datasets }, options: newChartOptions });
    }
}

async function plotSignalRelatedCharts(inputBufferForPlot = null, outputBufferForPlot = null) {
    const inputWaveCanvas = document.getElementById('inputWaveformChart');
    const outputWaveCanvas = document.getElementById('outputWaveformChart');
    if(!inputWaveCanvas || !outputWaveCanvas) return;

    const inputWaveCtx = inputWaveCanvas.getContext('2d', { willReadFrequently: true });
    const outputWaveCtx = outputWaveCanvas.getContext('2d', { willReadFrequently: true });

    const displayDuration = 0.05;
    const fs = audioContext ? audioContext.sampleRate : getSamplingFrequency();
    const signalAmplitude = getSignalAmplitude();
    const maxDisplaySamples = Math.floor(displayDuration * fs);

    let finalInputBuffer = inputBufferForPlot;
    let finalOutputBuffer = outputBufferForPlot;

    if (!finalInputBuffer) {
        const signalType = document.getElementById('signalType').value;
        finalInputBuffer = generateSignalBuffer(signalType, displayDuration, fs, signalAmplitude);
    }

    if (finalInputBuffer && !outputBufferForPlot) {
        const coeffs = getCoefficients();
        if (!audioContext) { 
            console.warn("AudioContext not available for offline rendering in plotSignalRelatedCharts.");
            finalOutputBuffer = null;
        } else {
            const offlineCtx = new OfflineAudioContext(1, finalInputBuffer.length, fs);
            const offlineSource = offlineCtx.createBufferSource();
            offlineSource.buffer = finalInputBuffer;
            const offlineIIRFilter = offlineCtx.createIIRFilter([coeffs.b0, coeffs.b1, coeffs.b2], [coeffs.a0, coeffs.a1, coeffs.a2]);
            offlineSource.connect(offlineIIRFilter);
            offlineIIRFilter.connect(offlineCtx.destination);
            offlineSource.start();
            try {
                const renderedBuffer = await offlineCtx.startRendering();
                const outputChannelData = renderedBuffer.getChannelData(0);
                for (let i = 0; i < outputChannelData.length; i++) {
                    outputChannelData[i] = Math.max(-1.0, Math.min(1.0, outputChannelData[i]));
                }
                finalOutputBuffer = audioContext.createBuffer(1, renderedBuffer.length, renderedBuffer.sampleRate);
                finalOutputBuffer.copyToChannel(outputChannelData, 0);
            } catch (e) {
                console.error("Error rendering output buffer for plotting: ", e);
                finalOutputBuffer = null;
            }
        }
    }

    const inputDataPoints = [];
    if (finalInputBuffer) {
        const inputChannelData = finalInputBuffer.getChannelData(0);
        const samplesToPlot = Math.min(inputChannelData.length, maxDisplaySamples);
        for (let i = 0; i < samplesToPlot; i++) {
            inputDataPoints.push({ x: (i / fs) * 1000, y: inputChannelData[i] });
        }
    }
    const outputDataPoints = [];
    if (finalOutputBuffer) {
        const outputChannelData = finalOutputBuffer.getChannelData(0);
        const samplesToPlot = Math.min(outputChannelData.length, maxDisplaySamples);
        for (let i = 0; i < samplesToPlot; i++) {
            outputDataPoints.push({ x: (i / fs) * 1000, y: outputChannelData[i] });
        }
    }

    if (inputWaveformChart) {
        inputWaveformChart.data.datasets[0].data = inputDataPoints;
        inputWaveformChart.options.scales.x.max = displayDuration * 1000;
        inputWaveformChart.options.scales.y.min = - (1);
        inputWaveformChart.options.scales.y.max = 1;
        inputWaveformChart.update('none');
    } else {
        inputWaveformChart = new Chart(inputWaveCtx, {
            type: 'line', data: { datasets: [{ label: '入力信号', data: inputDataPoints, borderColor: 'rgb(0, 123, 255)', tension: 0.1, pointRadius: 0 }] },
            options: { responsive: true, maintainAspectRatio: false, animation: false, scales: { x: { type: 'linear', title: { display: true, text: '時間 (ms)' }, min: 0, max: displayDuration * 1000 }, y: { title: { display: true, text: '振幅' }, min: -(1), max: 1 } } }
        });
    }
    if (outputWaveformChart) {
        outputWaveformChart.data.datasets[0].data = outputDataPoints;
        outputWaveformChart.options.scales.x.max = displayDuration * 1000;
        outputWaveformChart.options.scales.y.min = -(1);
        outputWaveformChart.options.scales.y.max = 1;
        outputWaveformChart.update('none');
    } else {
        outputWaveformChart = new Chart(outputWaveCtx, {
            type: 'line', data: { datasets: [{ label: '出力信号', data: outputDataPoints, borderColor: 'rgb(255, 100, 132)', tension: 0.1, pointRadius: 0 }] },
            options: { responsive: true, maintainAspectRatio: false, animation: false, scales: { x: { type: 'linear', title: { display: true, text: '時間 (ms)' }, min: 0, max: displayDuration * 1000 }, y: { title: { display: true, text: '振幅' }, min: -(1), max: 1 } } }
        });
    }

    inputFftChart = plotFftChart(inputFftChart, 'inputFftChart', finalInputBuffer, fs, '入力');
    outputFftChart = plotFftChart(outputFftChart, 'outputFftChart', finalOutputBuffer, fs, '出力');
}

async function updateAllPlots() {
    const coeffs = getCoefficients();
    const fs = getSamplingFrequency();

    const freqData = calculateFrequencyResponse(coeffs, fs);
    freqChart = plotFrequencyResponse(freqData, freqChart, 'freqResponseChart');

    const pzData = calculatePolesAndZeros(coeffs);
    plotZPlanePolesZeros(pzData);

    const impulseData = calculateImpulseResponse(coeffs);
    plotImpulseResponse(impulseData);

    if (!window.isPlotting3D) { 
        window.isPlotting3D = true; 
        try {
            const zPlane3DData = calculateZPlane3DData(coeffs);
            await plotZPlane3D(zPlane3DData); 
        } catch (error) {
            console.error("Error during 3D plot update:", error);
        } finally {
            window.isPlotting3D = false; 
        }
    }
    plotSignalRelatedCharts();
}
