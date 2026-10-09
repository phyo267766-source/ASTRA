# ASTRA Combined — React + Tailwind + Python Flask

This project combines the existing React/Tailwind UI and Flask SPHEREx search/compare backend. It also reconstructs a **safe exploratory region catalog** from the idea in the uploaded incomplete `main (2).py`. The uploaded Python fragment was **not executable as provided**; its simulated image fallback, invented Planet X classifications, and randomized historical light curves have deliberately NOT been included.

## Windows quick start

Requires Python 3.11+ and Node.js 20+ with npm. Internet is needed for dependencies and NASA requests.

### Terminal 1 (backend)
```powershell
cd backend
python -m pip install -r requirements.txt
python app.py
```

### Terminal 2 (frontend)
```powershell
cd frontend
npm install
npm run dev
```
Open http://127.0.0.1:5173 (not 5000) for the React interface. Vite proxies /api to Flask on port 5000.

## Layout
Single-page mission dashboard. No left navigation sidebar or right information sidebar. Search controls, sky viewer, observations, candidate catalog and recent activity are stacked on the same page.

## Features
- Observation search from NASA/IPAC IRSA; select two epochs and load real FITS cutouts
- WCS alignment, normalized difference overlay, and candidate-region list from connected flagged pixels
- No fabricated detections, magnitudes, or light curves; regions are **not verified astronomical candidates**
- A zero-result archive search may be caused by a lack of data or overly strict access_url filtering. Debug logs show raw vs eligible rows. This integration has not been verified against the live archive.

## Limitations
The FITS cutout API, access URL restrictions, and band compatibility must be validated against current IRSA responses. The exploratory 99th-percentile threshold guarantees some flagged pixels on most nonconstant comparisons; flagged regions do not establish moving bodies, photometric variability, or detection confidence. No historical light curve can be inferred from only two images.
