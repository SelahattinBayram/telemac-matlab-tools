;;; OUTLINETXT exports one closed 2-D polyline as an ordered X Y text file.
;;; Load with APPLOAD and run OUTLINETXT in AutoCAD for Mac or Windows.
;;; X and Y are written in World coordinates. The first point is repeated last.
;;; Circular polyline arcs are sampled with a user-selected maximum deviation.
;;; The drawing and the selected polyline are read only.

(defun otxt:value (code data fallback / item)
  (if (setq item (assoc code data)) (cdr item) fallback)
)

(defun otxt:stop (message)
  (princ (strcat "\n" message))
  (exit)
)

(defun otxt:vertices (entity data / kind records point bulge item vertex entry result)
  (setq kind (cdr (assoc 0 data)) result nil)
  (if (= kind "LWPOLYLINE")
    (progn
      (setq point nil bulge 0.0)
      (foreach item data
        (cond
          ((= (car item) 10)
            (if point (setq result (cons (list point bulge) result)))
            (setq point (cdr item) bulge 0.0))
          ((= (car item) 42) (setq bulge (cdr item)))
        )
      )
      (if point (setq result (cons (list point bulge) result)))
    )
    (progn
      (setq vertex (entnext entity))
      (while (and vertex (/= (cdr (assoc 0 (entget vertex))) "SEQEND"))
        (setq entry (entget vertex))
        (if (= (cdr (assoc 0 entry)) "VERTEX")
          (setq result
            (cons (list (cdr (assoc 10 entry)) (otxt:value 42 entry 0.0)) result))
        )
        (setq vertex (entnext vertex))
      )
    )
  )
  (reverse result)
)

(defun otxt:arc (a b bulge tolerance / dx dy chord factor cx cy radius sweep angle
                                         cosine step count index points phi)
  (setq dx (- (car b) (car a)) dy (- (cadr b) (cadr a))
        chord (sqrt (+ (* dx dx) (* dy dy))))
  (if (< chord 1e-10)
    nil
    (if (< (abs bulge) 1e-12)
      (list (list (car a) (cadr a)))
      (progn
        (setq factor (/ (- 1.0 (* bulge bulge)) (* 4.0 bulge))
              cx (- (/ (+ (car a) (car b)) 2.0) (* factor dy))
              cy (+ (/ (+ (cadr a) (cadr b)) 2.0) (* factor dx))
              radius (/ (* chord (+ 1.0 (* bulge bulge))) (* 4.0 (abs bulge)))
              sweep (* 4.0 (atan bulge))
              angle (atan (- (cadr a) cy) (- (car a) cx))
              cosine (max -1.0 (min 1.0 (- 1.0 (/ tolerance radius))))
              step (min (/ pi 2.0)
                     (* 2.0 (atan (sqrt (max 0.0 (- 1.0 (* cosine cosine)))) cosine))))
        (if (< step 1e-10)
          (otxt:stop "Arc tolerance is too small. Use a larger tolerance."))
        (setq count (max 1 (fix (+ 0.999999999999 (/ (abs sweep) step)))))
        (if (> count 100000)
          (otxt:stop "Arc requires too many points. Increase the tolerance."))
        (setq index 1 points (list (list (car a) (cadr a))))
        (while (< index count)
          (setq phi (+ angle (* sweep (/ (float index) count)))
                points (cons (list (+ cx (* radius (cos phi)))
                                   (+ cy (* radius (sin phi)))) points)
                index (1+ index))
        )
        (reverse points)
      )
    )
  )
)

(defun otxt:points (entity vertices elevation tolerance / remaining a b points point world result)
  (setq remaining vertices result nil)
  (while remaining
    (setq a (car remaining)
          b (if (cdr remaining) (cadr remaining) (car vertices))
          points (otxt:arc (car a) (car b) (cadr a) tolerance))
    (foreach point points
      (setq world (trans (list (car point) (cadr point) elevation) entity 0))
      (setq result (cons (list (car world) (cadr world)) result))
    )
    (setq remaining (cdr remaining))
  )
  (setq result (reverse result))
  (if (< (length result) 3)
    (otxt:stop "Outline must contain at least three distinct vertices."))
  (append result (list (car result)))
)

(defun c:OUTLINETXT (/ *error* fd selection entity data kind flags normal vertices
                      elevation tolerance hasArcs points filepath point)
  (defun *error* (message)
    (if fd (progn (close fd) (setq fd nil)))
    (if (and message (not (wcmatch (strcase message) "*BREAK*,*CANCEL*,*EXIT*")))
      (princ (strcat "\nOUTLINETXT error: " message)))
    (princ)
  )
  (setq selection (entsel "\nSelect the closed outline polyline to export as TXT: "))
  (if selection
    (progn
      (setq entity (car selection) data (entget entity)
            kind (cdr (assoc 0 data)) flags (otxt:value 70 data 0))
      (if (not (member kind '("LWPOLYLINE" "POLYLINE")))
        (otxt:stop "Select a closed 2D polyline. Draw it with PLINE and use Close."))
      (if (= (logand flags 1) 0)
        (otxt:stop "Polyline is open. Use PEDIT > Close and try again."))
      (if (and (= kind "POLYLINE") (/= (logand flags 94) 0))
        (otxt:stop "3D, mesh or spline-fit polylines are not supported. Use a planar 2D PLINE."))
      (setq normal (otxt:value 210 data '(0.0 0.0 1.0)))
      (if (or (> (abs (car normal)) 1e-10) (> (abs (cadr normal)) 1e-10))
        (otxt:stop "Polyline must be parallel to the XY plane."))
      (setq elevation
        (if (= kind "LWPOLYLINE")
          (otxt:value 38 data 0.0)
          (caddr (otxt:value 10 data '(0.0 0.0 0.0)))))
      (setq vertices (otxt:vertices entity data) hasArcs nil tolerance 0.5)
      (if (< (length vertices) 3)
        (otxt:stop "Outline must contain at least three vertices."))
      (foreach point vertices
        (if (> (abs (cadr point)) 1e-12) (setq hasArcs T)))
      (if hasArcs
        (progn
          (initget 6)
          (setq tolerance (getreal "\nMaximum arc sampling deviation [drawing units] <0.5>: "))
          (if (null tolerance) (setq tolerance 0.5))
        )
      )
      (setq points (otxt:points entity vertices elevation tolerance))
      (setq filepath
        (getfiled "Save outline TXT"
          (strcat (getvar "DWGPREFIX") "outline_autocad.txt") "txt" 1))
      (if filepath
        (progn
          (setq fd (open filepath "w"))
          (if (null fd) (otxt:stop "Cannot open TXT file. Check the output location."))
          (foreach point points
            (write-line (strcat (rtos (car point) 2 10) (chr 9)
                                (rtos (cadr point) 2 10)) fd))
          (close fd)
          (setq fd nil)
          (princ (strcat "\nOutline saved: " filepath))
          (princ (strcat "\n" (itoa (length points))
                         " XY rows. The first point is repeated at the end to close the boundary."))
          (princ "\nSelect this file as the outer boundary in MATLAB telemac_preprocess.")
        )
        (princ "\nSave cancelled.")
      )
    )
    (princ "\nSelection cancelled.")
  )
  (princ)
)

(princ "\nOutline exporter loaded. Command: OUTLINETXT")
(princ)
