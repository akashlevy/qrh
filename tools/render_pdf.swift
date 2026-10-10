// Render every page of a PDF to PNG with PDFKit: swift render_pdf.swift in.pdf outdir [dpi]
import Foundation
import PDFKit
import AppKit
let args = CommandLine.arguments
guard args.count >= 3, let doc = PDFDocument(url: URL(fileURLWithPath: args[1])) else {
  fputs("usage: render_pdf.swift in.pdf outdir [dpi]\n", stderr); exit(2)
}
let dpi = args.count > 3 ? Double(args[3])! : 110
let out = URL(fileURLWithPath: args[2])
try? FileManager.default.createDirectory(at: out, withIntermediateDirectories: true)
for i in 0..<doc.pageCount {
  let page = doc.page(at: i)!
  let box = page.bounds(for: .mediaBox)
  let w = Int(box.width * dpi / 72), h = Int(box.height * dpi / 72)
  let rep = NSBitmapImageRep(bitmapDataPlanes: nil, pixelsWide: w, pixelsHigh: h, bitsPerSample: 8,
    samplesPerPixel: 4, hasAlpha: true, isPlanar: false, colorSpaceName: .deviceRGB, bytesPerRow: 0, bitsPerPixel: 0)!
  let ctx = NSGraphicsContext(bitmapImageRep: rep)!.cgContext
  ctx.setFillColor(CGColor(red: 1, green: 1, blue: 1, alpha: 1)); ctx.fill(CGRect(x: 0, y: 0, width: w, height: h))
  ctx.scaleBy(x: CGFloat(dpi / 72), y: CGFloat(dpi / 72))
  page.draw(with: .mediaBox, to: ctx)
  let png = rep.representation(using: .png, properties: [:])!
  try! png.write(to: out.appendingPathComponent(String(format: "page-%02d.png", i + 1)))
}
print("pages: \(doc.pageCount)")
